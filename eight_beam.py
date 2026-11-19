"""Same experiment protocol, independently declared eight-beam geometry."""
from pathlib import Path
import study
if __name__=='__main__':
    study.GEOMETRY='oblique8'
    study.DIM=8
    study.OUT=Path('results/eight_beam')
    study.OUT.mkdir(exist_ok=True)
    study.main()
