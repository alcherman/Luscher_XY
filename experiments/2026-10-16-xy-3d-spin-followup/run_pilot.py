"""Finite launchd pilot supervisor; no new long runs without measured review."""
from pathlib import Path
import sys,json,subprocess,fcntl,os,time,traceback
import campaign as c
P=c.OUT

def main():
    with (P/'pilot-manager.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);awake=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(os.getpid())])
        state=dict(pid=os.getpid(),started=c.frozen.now(),stages={})
        def publish(value,**extra):state.update(state=value,updated=c.frozen.now(),**extra);c.atomic(P/'pilot-status.json',state)
        try:
            config=c.read(P/'control.json');config['launch_new']=True;c.atomic(P/'control.json',config)
            for tag,workers in [('smoke',1),('benchmark',1),('pilot',8)]:
                publish('running_'+tag)
                with (P/'logs'/f'coordinator-{tag}.log').open('a') as stream:
                    subprocess.run([sys.executable,'-u',str(P/'campaign.py'),'--coordinate','--manifest',str(P/'batches'/f'{tag}.json'),'--workers',str(workers)],cwd=c.ROOT,stdout=stream,stderr=subprocess.STDOUT,check=True)
                batch=c.read(P/'batches'/f'{tag}.json');jobs=[c.read(P/'jobs'/f"{p['id']}.json") for p in batch['points']]
                if any(j['state']!='complete' or not j.get('audit',{}).get('passed') or (j['point']['online'] and not j['audit'].get('online_spin',{}).get('passed')) for j in jobs):raise RuntimeError('Closure audit failed:'+tag)
                state['stages'][tag]=dict(jobs=len(jobs),generation_worker_seconds=sum(j['generation_seconds'] for j in jobs),audit_worker_seconds=sum(j['audit_seconds'] for j in jobs),full_bytes=sum(j['audit']['bytes'] for j in jobs));publish('closed_'+tag)
                if tag=='benchmark':
                    import h5py,numpy as np
                    with h5py.File(P/'chains'/f"{jobs[0]['point']['id']}.h5") as a,h5py.File(P/'chains'/f"{jobs[1]['point']['id']}.h5") as b:
                        for key in a:
                            if not key.endswith('_seconds'):np.testing.assert_array_equal(a[key][:],b[key][:],err_msg=key)
                    state['benchmark_physical_rng_histories_exact']=True
            config=c.read(P/'control.json');config['launch_new']=False;c.atomic(P/'control.json',config);publish('pilot_audited_analysing')
            with (P/'logs/pilot-analysis.log').open('a') as stream:subprocess.run([sys.executable,'-u',str(P/'analyze_pilot.py')],cwd=c.ROOT,stdout=stream,stderr=subprocess.STDOUT,check=True)
            publish('pilot_complete_pending_main_design_review')
        except BaseException:
            config=c.read(P/'control.json');config['launch_new']=False;c.atomic(P/'control.json',config);publish('failed_needs_attention',error=traceback.format_exc());raise
        finally:
            if awake.poll() is None:awake.terminate()
if __name__=='__main__':main()
