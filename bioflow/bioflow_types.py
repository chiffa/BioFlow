from typing import NewType
from collections.abc import Callable
from collections.abc import Mapping, Sequence
from collections.abc import Iterable
from typing import Any, Union, Tuple, List, Dict


Optionally_Weighted_Internal_IDs = NewType('Optionally_Weighted_Internal_IDs',
                                           Union[List[int], List[Tuple[int, float]]])

# TODO:
#   - annotome sample entry
#   - interactome sample entry