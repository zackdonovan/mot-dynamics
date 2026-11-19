import json
from pathlib import Path
from spatial_motion import run
r,*_=run('oblique8',[.006,.003,.002],T=8000.,rtol=2e-10)
Path('results/recoil_validity/long_tolerance_check.json').write_text(json.dumps(r,indent=2))
