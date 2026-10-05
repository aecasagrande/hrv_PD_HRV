"""
Run the HRV pipeline (`main`) on the synthetic recordings in ./mock_data
without editing it: loads the script as a module, points BASE / OUTPUT_DIR at
./mock_data and ./output, sets BATCH_MODE = True and DRY_RUN = False, and runs
the batch over the whole PARTICIPANTS registry (entries without mock data are
reported as 'no trials'). Then checks the workbooks against the ground truth.

    python generate_mock_data.py
    python run_simulation.py
"""
import importlib.machinery
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

if __name__ == '__main__':
    os.chdir(HERE)
    loader = importlib.machinery.SourceFileLoader('hrv_pipeline', os.path.join(HERE, 'main'))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    hrv = importlib.util.module_from_spec(spec)
    loader.exec_module(hrv)
    hrv.BASE, hrv.OUTPUT_DIR = './mock_data', './output'
    hrv.BATCH_MODE, hrv.DRY_RUN = True, False
    hrv.main()

    import validate_mock_results
    sys.exit(0 if validate_mock_results.main() else 1)
