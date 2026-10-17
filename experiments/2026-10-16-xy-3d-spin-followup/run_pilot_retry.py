"""Finite retry of pilot after pre-generation target guard failure."""
import os,sys,fcntl,subprocess,traceback
import campaign as c
P=c.OUT

def main():
    with (P/'pilot-manager.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        awake=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(os.getpid())])
        state=dict(pid=os.getpid(),started=c.frozen.now(),manifest='pilot-retry-v2.json',prior_status='attempt-history/pilot-status-warmup-guard.json')
        def publish(value,**extra):
            state.update(state=value,updated=c.frozen.now(),**extra);c.atomic(P/'pilot-status.json',state)
        try:
            config=c.read(P/'control.json');config['launch_new']=True;c.atomic(P/'control.json',config)
            publish('running_pilot_retry')
            with (P/'logs/coordinator-pilot-retry-v2.log').open('a') as stream:
                subprocess.run([sys.executable,'-u',str(P/'campaign.py'),'--coordinate','--manifest',str(P/'batches/pilot-retry-v2.json'),'--workers','8'],cwd=c.ROOT,stdout=stream,stderr=subprocess.STDOUT,check=True)
            batch=c.read(P/'batches/pilot-retry-v2.json');jobs=[c.read(P/'jobs'/f"{p['id']}.json") for p in batch['points']]
            if any(j['state']!='complete' or not j.get('audit',{}).get('passed') or not j['audit'].get('online_spin',{}).get('passed') for j in jobs):raise RuntimeError('Pilot closure audit failed')
            config=c.read(P/'control.json');config['launch_new']=False;c.atomic(P/'control.json',config)
            publish('pilot_audited_analysing',audited_jobs=len(jobs))
            with (P/'logs/pilot-analysis-retry-v2.log').open('a') as stream:
                subprocess.run([sys.executable,'-u',str(P/'analyze_pilot.py')],cwd=c.ROOT,stdout=stream,stderr=subprocess.STDOUT,check=True)
            publish('pilot_complete_pending_main_design_review')
        except BaseException:
            config=c.read(P/'control.json');config['launch_new']=False;c.atomic(P/'control.json',config)
            publish('failed_needs_attention',error=traceback.format_exc());raise
        finally:
            if awake.poll() is None:awake.terminate()
if __name__=='__main__':main()
