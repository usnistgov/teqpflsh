import subprocess, os, shutil, re, timeit

here = os.path.dirname(__file__)
import teqpflsh
import argparse 

def run(doxygen=True, notebooks=False):
    
    if doxygen:
        if os.path.exists(here+'/_static/'):
            shutil.rmtree(here+'/_static/')
        os.makedirs(here+'/_static')

        oldDoxyfile = here+'/../Doxyfile'
        newDoxyfile = here+'/../../Doxyfile.injected'
        def repl(matchobj):
            return matchobj.group(1) + teqpflsh.__version__
        newcontents = re.sub(r'(PROJECT_NUMBER\s+=\s+)(.+)', repl, open(oldDoxyfile).read())
        with open(newDoxyfile, 'w') as fp:
            fp.write(newcontents)

        subprocess.check_call('doxygen Doxyfile.injected', cwd=here+'/../..', shell=True)
        os.remove(newDoxyfile)
        
    if notebooks:
        # Execute all the notebooks to populate all cells
        for path, dirs, files in os.walk(here):
            for file in files:
                if file.endswith('.ipynb') and '.ipynb_checkpoints' not in path:
                    print(f'Converting {file} in {path}')
                    tic = timeit.default_timer()
                    subprocess.check_output(f'jupyter nbconvert --allow-errors --to notebook --output "{file}" --execute "{file}"', shell=True, cwd=path)
                    toc = timeit.default_timer()
                    print(f'Elapsed for conversion: {toc-tic} seconds')

if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description='Pre-runner.')
    parser.add_argument('--no-doxygen', action='store_true', help="Run doxygen")
    parser.add_argument('--no-notebooks', action='store_true', help="Execute the Jupyter notebooks")
    args = parser.parse_args()
    
    print('running sphinx_pre_run.py...')
    run(doxygen=not args.no_doxygen, notebooks=not args.no_notebooks)
