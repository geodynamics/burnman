import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
import glob, os

stream = os.popen('jupyter kernelspec list')
output = stream.read()

try:
    highest_kernel_name = output.strip().split()[-2]
except IndexError as e:
    print('Jupyter kernel name not found.')
    exit()

for notebook_filename in sorted(glob.glob("../tutorial/*.ipynb")):
    print(f'Running {notebook_filename}')

    with open(notebook_filename) as f:
        nb = nbformat.read(f, as_version=4)
        ep = ExecutePreprocessor(timeout=600, kernel_name=highest_kernel_name)
        ep.preprocess(nb, {'metadata': {'path': ''}})

    print('Text outputs:')
    for cell in nb['cells']:
        try:
            for output in cell['outputs']:
                if output['name'] == 'stdout':
                    print(output['text'])
        except KeyError:
            pass

    print(f'Successfully finished running {notebook_filename}')
