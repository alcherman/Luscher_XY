"""Read-only long-distance spin measurements on uniform retained checkpoints.

No simulation calls; no adaptive/dense-window sampling in central estimates.
Every spatial average is one MC-time observation, not independent sites.
"""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import argparse,json,time,hashlib,fcntl,shutil
import h5py,numpy as np
OUT=Path(__file__).resolve().parent
REPO=OUT.parents[1]
CAM=OUT.parent/'2026-10-13-xy-infinite-volume'
TRANSFER=OUT.parent/'2026-10-16-mc-data-hdd-archive'
ARCHIVE=Path('/Volumes/MC_Data_Backup/Luscher_XY/archive-2026-10-16')
PREFIX='raw/Luscher_XY/2026-10-13-xy-infinite-volume/chains/'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atomic(p,a):
    tmp=p.with_suffix(p.suffix+'.tmp')
    with tmp.open('w') as f:json.dump(a,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(tmp,p)
def signature(path):
    s=Path(path).stat();return [s.st_size,s.st_mtime_ns,s.st_ino,s.st_dev]

def spin_measure(phi):
    raw=np.asarray(phi)
    if raw.ndim==4 and raw.shape[0]==1:raw=raw[0]
    if raw.ndim!=3 or len(set(raw.shape))!=1 or raw.shape[0]%4:raise ValueError('Need a cubic box with L divisible by4')
    if not np.isfinite(raw).all():raise ValueError('Missing/nonfinite retained field')
    L=raw.shape[0];V=L**3;spin=np.exp(1j*raw);axis=np.empty((2,3),complex)
    for k,r in enumerate([L//4,L//2]):
        for d in range(3):axis[k,d]=np.mean(np.conj(spin)*np.roll(spin,-r,axis=d))
    mag=spin.mean();m2=float(abs(mag)**2);f=[]
    phase=np.exp(-2j*np.pi*np.arange(L)/L)
    for d in range(3):
        profile=spin.sum(axis=tuple(a for a in range(3) if a!=d));f.append(float(abs(profile@phase)**2/V))
    return dict(axis_real=axis.real,axis_imag=axis.imag,m2=m2,S0=V*m2,F=np.array(f))

def retained_checkpoint_rows(h,steps,stride=100):
    if len(h['phi'])!=steps+1 or len(h['sweep'])!=steps+1 or int(h.attrs['committed_rows'])!=steps+1:raise ValueError('Incomplete/mismatched chain')
    rows=np.arange(0,steps+1,stride,dtype=np.int64)
    if steps%stride:raise ValueError('Final cycle not on prescribed cadence')
    if 'retained_indices' in h['phi'].attrs and not set(rows.tolist())<=set(map(int,h['phi'].attrs['retained_indices'])):raise ValueError('A prescribed regular checkpoint is absent')
    if not np.array_equal(h['sweep'][rows],rows):raise ValueError('Saved cycle numbers do not match field row indices')
    return rows

def freeze_plan(path):
    source=OUT/'joint-fss-existing-20261016-full/sources.json';groups=OUT/'joint-fss-existing-20261016-full/input-groups.json'
    a=json.loads(source.read_text());jobs={v['id']:v for v in a['source_jobs']};points=[]
    verified=json.loads((TRANSFER/'verified.json').read_text());wanted={i for g in json.loads(groups.read_text()) if g['L'] in [24,32,48,64] for i in g['ids']}
    checks={}
    with (TRANSFER/'checksums.jsonl').open() as f:
        for line in f:
            r=json.loads(line)
            if r['relative'].startswith(PREFIX) and Path(r['relative']).stem in wanted:checks[Path(r['relative']).stem]=r
    for g in json.loads(groups.read_text()):
        if g['L'] not in [24,32,48,64]:continue
        for i in g['ids']:
            p=jobs[i]['point'];relative=PREFIX+i+'.h5';src=ARCHIVE/relative
            if p['D']!=3 or p['beta']!=0:raise ValueError('Wrong model coordinate')
            if i not in checks or relative not in verified:raise ValueError('Missing completed transfer verification')
            if signature(src)!=verified[relative]:raise ValueError('Archive identity differs from verified transfer')
            points.append(dict(point=p,group_q=g['q'],actual_q=p['epsilon']/np.pi,source_path=str(src),source_signature=signature(src),archive_verified_sha256=checks[i]['sha256'],archive_sha_scope='Inherited from complete read-back SHA256 verification during transfer; whole file is not rehashed by this analysis. Current inode/size/mtime pinned; selected field bytes hashed and stored moments independently recomputed.'))
    plan=dict(created=time.time(),scope='Existing archived fields only, read-only. No new MC, compaction or environment changes.',points=points,pilot_ids=[v['point']['id'] for v in points if abs(v['group_q']-.80075)<1e-12],sizes=[24,32,48,64],separations=['L/4','L/2'],directions='All3 positive axes; retain individual complex-axis correlations and average REAL parts.',cadence_cycles=100,analysis_cut_cycles=1000,source_sha256=sha(source),source_groups_sha256=sha(groups),transfer_verified_sha256=sha(TRANSFER/'verified.json'),regular_sampling='Every100cycles, including warmup in saved derived history. Central estimates use cycles>1000; dense diagnostic/adaptively retained fields not pooled.',statistical_plan='Use matching-checkpoint m2,S0,F,winding2 and G histories. Blocks expressed in original cycle units: max dense monitored tau and100*G/checkpoint tau; do not regard checkpoints as independent. Log applied only to ensemble means/draws; preserve invalid tails. Pilot before expansion; no fixed-q slope interpreted as asymptotic Delta.')
    if len(points)!=64 or len(plan['pilot_ids'])!=8:raise ValueError('Unexpected source coverage')
    if path.exists():raise FileExistsError(path)
    atomic(path,plan)


def extract_one(entry,folder):
    p=entry['point'];i=p['id'];src=Path(entry['source_path']);before=signature(src)
    if before!=entry['source_signature']:raise ValueError('Archive changed since plan: '+i)
    lock=(CAM/'jobs'/(i+'.lock')).open('a+')
    fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
    start=time.monotonic();series=[];axes=[];imag=[];field_hash=hashlib.sha256();maxerr=np.zeros(3);policy={};params={}
    try:
        with h5py.File(src,'r') as h:
            handle=h.id.get_vfd_handle();fd=os.fstat(handle)
            if [fd.st_size,fd.st_mtime_ns,fd.st_ino,fd.st_dev]!=before:raise ValueError('Path/open-file identity mismatch')
            params=json.loads(h.attrs['parameters'])
            if params['D']!=3 or params['N']!=p['N'] or params['kappa']!=0 or abs(params['epsilon']-p['epsilon'])>1e-14:raise ValueError('Archived parameters mismatch')
            rows=retained_checkpoint_rows(h,p['steps']);stored_m=h['magnetization_squared'][rows];stored_z=h['structure_zero'][rows];stored_f=h['structure_min'][rows]
            wind=h['winding'][rows].astype(float);w2=(wind**2).mean(axis=(1,2));constant_winding_axes=[bool(np.all(wind[:,a,:]==wind[0,a,:])) for a in range(3)]
            if 'retention_record' in h.attrs:
                retention=json.loads(h.attrs['retention_record']);policy=retention.get('policy',{})
            for k,row in enumerate(rows):
                raw=h['phi'][int(row)];field_hash.update(np.array([row],dtype='<i8').tobytes());field_hash.update(np.ascontiguousarray(raw).tobytes());m=spin_measure(raw)
                np.testing.assert_allclose(m['m2'],stored_m[k],rtol=1e-10,atol=1e-12)
                np.testing.assert_allclose(m['S0'],stored_z[k],rtol=1e-10,atol=1e-8)
                np.testing.assert_allclose(m['F'],stored_f[k],rtol=1e-10,atol=1e-8)
                if abs(m['axis_imag'][1]).max()>1e-12:raise ValueError('Half-box sine reversal identity failed')
                maxerr=np.maximum(maxerr,[abs(m['m2']-stored_m[k]),abs(m['S0']-stored_z[k]),float(abs(m['F']-stored_f[k]).max())]);axes.append(m['axis_real']);imag.append(m['axis_imag']);series.append([*m['axis_real'].mean(axis=1),m['m2'],m['S0'],float(m['F'].mean()),w2[k]])
            # Dense scalar moments are cheap; record them to diagnose time sampling,
            # without substituting their estimates into the matched-checkpoint bootstrap.
            dense=dict(sweep=h['sweep'][:],m2=h['magnetization_squared'][:],S0=h['structure_zero'][:],F=h['structure_min'][:].mean(axis=1))
            if signature(src)!=before:raise ValueError('Archive changed during extraction')
    finally:lock.close()
    npz=folder/(i+'.npz')
    if npz.exists() or (folder/(i+'.json')).exists():raise FileExistsError('Refuse to overwrite derived measurements')
    np.savez_compressed(npz,sweep=rows,series=np.array(series),axis_real=np.array(axes),axis_imag=np.array(imag),**{'dense_'+k:v for k,v in dense.items()})
    meta=dict(id=i,point=p,group_q=entry['group_q'],actual_q=entry['actual_q'],source_path=str(src),source_signature=before,archive_verified_sha256=entry['archive_verified_sha256'],archive_sha_scope=entry['archive_sha_scope'],selected_fields_sha256=field_hash.hexdigest(),selected_cycle_count=len(rows),selected_postwarmup_count=int(np.sum(rows>1000)),selected_cycles=rows.tolist(),configuration_observation='Each G value is a spatial/axis mean of ONE configuration, not independent sites.',columns=['G_quarter','G_half','m2','S0','F','winding2'],max_abs_recompute_errors=dict(zip(['m2','S0','F_axis'],maxerr.tolist())),half_box_sine_identity_passed=True,original_retention_tau_estimates_cycles=policy.get('tau_estimates',{}),constant_winding_axes_all_checkpoints=constant_winding_axes,run_parameters=params,runtime_seconds=time.monotonic()-start,field_bytes_read=int(len(rows)*8*p['N']**3),derived_data_sha256=sha(npz))
    atomic(folder/(i+'.json'),meta)
    print(i,'fields',len(rows),'seconds',round(meta['runtime_seconds'],3),'moment_errors',maxerr.tolist(),flush=True)
    return meta


def main(args):
    planpath=Path(args.plan)
    if args.freeze:freeze_plan(planpath);return
    plan=json.loads(planpath.read_text());folder=Path(args.output);folder.mkdir(exist_ok=False)
    shutil.copyfile(__file__,folder/Path(__file__).name);atomic(folder/'design.json',dict(plan_path=str(planpath),plan_sha256=sha(planpath),phase=args.phase,implementation_sha256=sha(__file__)))
    points=[p for p in plan['points'] if (p['point']['id'] in plan['pilot_ids'])==(args.phase=='pilot')] if args.phase in ['pilot','remaining'] else plan['points']
    for n,entry in enumerate(points):
        extract_one(entry,folder);atomic(folder/'status.json',dict(state='extracting',completed=n+1,total=len(points),latest=entry['point']['id']))
    atomic(folder/'status.json',dict(state='measurements_complete_pending_statistical_audit',completed=len(points),total=len(points)))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--plan',required=True);ap.add_argument('--freeze',action='store_true');ap.add_argument('--phase',choices=['pilot','remaining','all'],default='pilot');ap.add_argument('--output');main(ap.parse_args())
