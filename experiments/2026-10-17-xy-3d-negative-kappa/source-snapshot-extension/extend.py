"""Second bounded finite-volume stage: move missed brackets downward in epsilon."""
import scout as s
import concurrent.futures,subprocess,sys,os,fcntl,traceback,json
P=s.P
NAME='continuation-bracket-extension.json'
def points():
    ps=[]
    for ki,k,qs in [(2,-.25,[.600,.650,.690]),(3,-.5,[.520,.570,.620])]:
        for q in qs:
            for n in [8,12]:
                for r in [2,3]:
                    seed=int(s.np.random.SeedSequence([2026101702,ki,round(q*1000),n,r]).generate_state(1,dtype=s.np.uint64)[0])
                    ps.append(dict(id=f'negative3_ext_k{ki}_q{round(q*1000):04d}_N{n:03d}_r{r}',D=3,N=n,beta=k,kappa=k,epsilon_over_pi=q,epsilon=s.np.pi*q,seed=seed,initial_winding=[r%2,0,0],steps=4000,warmup=1000,update='metropolis'))
    return ps

def main():
    with (P/'manager.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        review=s.read(P/'extension-review.json')
        if not review.get('approved') or review['source_hashes']!=s.hashes() or review['extension_sha256']!=s.legacy.sha256(__file__):raise RuntimeError('Frozen extension review required')
        ps=points();manifest=P/'batches'/NAME
        if not manifest.exists():
            guard=s.storage.acquire()
            try:
                copies=0
                for folder in [s.OLD,s.ROOT/'experiments/2026-10-16-xy-3d-spin-followup',P]:
                    s.storage.OUT=folder;copies+=s.storage.reserved_remaining()
                    for jf in (folder/'jobs').glob('*.json'):
                        if s.read(jf).get('state') in ['running','starting']:raise RuntimeError('Known active campaign requires reconciled storage before extension')
                free=s.shutil.disk_usage(P/'chains').free;ceil=sum(s.storage.estimated_full_bytes(p) for p in ps)
                if free-ceil-copies<25*2**30:raise RuntimeError('Full remaining budget fails')
                s.atomic(manifest,dict(points=ps,source_hashes=s.hashes(),extension_sha256=s.legacy.sha256(__file__),full_ceiling_bytes=ceil,physical_free_bytes=free,created=s.legacy.now(),scope='Lowerepsilon missed-bracket extension only, not independent infinite-volume curve'))
            finally:guard.close()
        elif s.read(manifest)['points']!=ps or s.read(manifest)['source_hashes']!=s.hashes():raise RuntimeError('Frozen extension differs')
        awake=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(os.getpid())])
        try:
            def launch(p):
                with (P/'logs'/f"{p['id']}.log").open('a') as out:subprocess.run([sys.executable,'-u',str(P/'scout.py'),'worker',p['id'],NAME],cwd=s.ROOT,stdout=out,stderr=subprocess.STDOUT,check=True)
                return p['id']
            s.atomic(P/'status.json',dict(state='running_lower_epsilon_extension',total=24,updated=s.legacy.now()))
            with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
                for i,identifier in enumerate(pool.map(launch,ps),1):s.atomic(P/'status.json',dict(state='running_lower_epsilon_extension',complete=i,total=24,last=identifier,updated=s.legacy.now()))
            with (P/'logs/analysis-v3.log').open('a') as out:subprocess.run([sys.executable,'-u',str(P/'analyse.py')],cwd=s.ROOT,stdout=out,stderr=subprocess.STDOUT,check=True)
            s.atomic(P/'status.json',dict(state='complete_preliminary_finite_volume',complete=72,updated=s.legacy.now()))
        except BaseException:
            s.atomic(P/'status.json',dict(state='failed_needs_attention',error=traceback.format_exc(),updated=s.legacy.now()));raise
        finally:
            if awake.poll() is None:awake.terminate()
if __name__=='__main__':main()
