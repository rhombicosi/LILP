"""Energy parameter set used by the loop classes, chosen once at import time:

    LILP_PARAMS=tnn         Turner 2004 (NNDB), default    -> LILP/tnn_parameters/
    LILP_PARAMS=andronescu  Andronescu et al. 2010, BL*    -> LILP/andronescu_parameters/

Both packages define the same tables under the same names; each module of this
package re-exports the matching module of the selected set. Loop classes import
from here (LILP.parameters.<module>), never from a set directly.
"""
import importlib
import os

PARAMETER_SETS = {
    'tnn': 'LILP.tnn_parameters',
    'andronescu': 'LILP.andronescu_parameters',
}

NAME = os.environ.get('LILP_PARAMS', 'tnn')
if NAME not in PARAMETER_SETS:
    raise ValueError(f'LILP_PARAMS={NAME!r}, expected one of {sorted(PARAMETER_SETS)}')


def export(module_name: str, namespace: dict) -> None:
    # copies the public names of <selected set>.<module_name> into namespace
    module = importlib.import_module(f'{PARAMETER_SETS[NAME]}.{module_name}')
    namespace.update({k: v for k, v in vars(module).items() if not k.startswith('_')})
