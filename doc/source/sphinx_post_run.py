import subprocess, os
                
def remove_cell_outputs():
    # Remove all outputs from all cells (to avoid large commits from the figures in cell outputs)
    for path, dirs, files in os.walk('.'):
        for file in files:
            if file.endswith('.ipynb') and '.ipynb_checkpoints' not in path:
                # subprocess.check_output(f'git checkout "{file}"', shell=True, cwd=path)
                subprocess.check_output(f'jupyter nbconvert --clear-output --inplace "{file}"', shell=True, cwd=path)

if __name__ == '__main__':
    print('running sphinx_post_run.py...')
    remove_cell_outputs()