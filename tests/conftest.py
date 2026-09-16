# icepack2_tools is not an installed package, it is imported from the repo
# root; this is what lets a bare `pytest tests/` from the repo root find it.
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
