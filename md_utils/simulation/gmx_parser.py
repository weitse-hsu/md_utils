import os
import re
import six
import warnings
from collections import OrderedDict


def parse_ndx(ndx_file):
    """
    Parse a GROMACS index (.ndx) file and return a dictionary of groups.

    Parameters
    ----------
    ndx_file : str
        Path to the GROMACS index file.

    Returns
    -------
    groups : dict
        An ordered dictionary where keys are group names and values are lists of atom indices.
    group_str : str
        A string representation of the groups for easy viewing.
    """
    groups = OrderedDict()
    current_group = None
    with open(ndx_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('[') and line.endswith(']'):
                current_group = line[1:-1].strip()
                groups[current_group] = []
            elif current_group is not None:
                indices = line.split()
                groups[current_group].extend(int(idx) for idx in indices)

    name_w = max((len(name) for name in groups), default=0)
    label_w = name_w + 1
    count_w = len(str(max((len(atoms) for atoms in groups.values()), default=0)))
    idx_w = len(str(len(groups) - 1))

    group_str = ''
    for i, (name, atoms) in enumerate(groups.items()):
        label = f"{name}:"
        group_str += f"  ({i}) {'':<{idx_w - len(str(i))}}{label:<{label_w}} {len(atoms):>{count_w}} atoms\n"

    return groups, group_str


class ParseError(Exception):
    """Error raised during parsing a file."""


class MDP(OrderedDict):
    """
    A class that represents a GROMACS MDP file. Note that an MDP instance is an ordered dictionary,
    with the i-th key corresponding to the i-th line in the MDP file. Comments and blank lines are
    also preserved, e.g., with keys 'C0001' and 'B0001', respectively. The value corresponding to a
    'C' key is the comment itself, while the value corresponding to a 'B' key is an empty string.
    Comments after a parameter on the same line are discarded. Leading and trailing spaces
    are always stripped.

    Parameters
    ----------
    input_mdp : str, Optional
        The path of the input MDP file. The default is None.
    **kwargs : Optional
        Additional keyword arguments to be passed to add additional key-value pairs to the MDP instance.
        Note that no sanity checks will be performed for the key-value pairs passed in this way. This
        also does not work for keys that are not legal python variable names, such as anything that includes
        a minus '-' sign or starts with a number.

    Attributes
    ----------
    COMMENT : :code:`re.Pattern` object
        A compiled regular expression pattern for comments in MDP files.
    PARAMETER : :code:`re.Pattern` object
        A compiled regular expression pattern for parameters in MDP files.
    input_mdp : str
        The real path of the input MDP file returned by :code:`os.path.realpath(input_mdp)`,
        which resolves any symbolic links in the path.

    Example
    -------
    >>> from ensemble_md.utils import gmx_parser
    >>> gmx_parser.MDP("em.mdp")
    MDP([('C0001', 'em.mdp - used as input into grompp to generate em.tpr'),
         ('integrator', 'steep'), ('nsteps', 500000), ...])
    """
    # Below are some class variables accessible to all functions.
    COMMENT = re.compile(r"""\s*;\s*(?P<value>.*)""")
    PARAMETER = re.compile(
        r"""\s*(?P<parameter>[^=]+?)\s*=\s*(?P<value>[^;]*)(?P<comment>\s*;.*)?""",
        re.VERBOSE,
    )

    def __init__(self, input_mdp=None, **kwargs):
        super(MDP, self).__init__(**kwargs)  # can use kwargs to set dict! (but no sanity checks!)
        if input_mdp is not None:
            self.input_mdp = os.path.realpath(input_mdp)
            self.read()

    def _convert_to_numeric(self, s):
        """
        Converts the input to a numerical type when possible. This internal function is used for the MDP parser.

        Parameters
        ----------
        s : any
            The input value to be converted to a numerical type if possible. The data type of :code:`s` is
            usually :code:`str` but can be any. However, if :code:`s` is not a string, it will be returned as is.

        Returns
        -------
        numerical : any
            The converted numerical value. If :code:`s` can be converted to a single numerical value,
            that value is returned as an :code:`int` or :code:`float`. If :code:`s` can be converted to
            multiple numerical values, a list containing those values is returned.
            If :code:`s` cannot be converted to a numerical value, :code:`s` is returned as is.
        """
        if type(s) is not str:
            return s
        for converter in int, float, str:  # try them in increasing order of lenience
            try:
                s = [converter(i) for i in s.split()]
                if len(s) == 1:
                    return s[0]
                else:
                    return s
            except (ValueError, AttributeError):
                pass

    def read(self):
        """
        Reads and parses the input MDP file.
        """
        def BLANK(i):
            return f"B{i:04d}"

        def COMMENT(i):
            return f"C{i:04d}"

        data = OrderedDict()
        iblank = icomment = 0
        with open(self.input_mdp) as mdp:
            for line in mdp:
                line = line.strip()
                if len(line) == 0:
                    iblank += 1
                    data[BLANK(iblank)] = ""
                    continue
                m = self.COMMENT.match(line)
                if m:
                    icomment += 1
                    data[COMMENT(icomment)] = m.group("value")
                    continue

                m = self.PARAMETER.match(line)
                if m:
                    parameter = m.group("parameter")
                    value = self._convert_to_numeric(m.group("value"))
                    data[parameter] = value
                else:
                    err_msg = f"{os.path.basename(self.input_mdp)!r}: unknown line in mdp file, {line!r}"
                    raise ParseError(err_msg)

        super(MDP, self).update(data)

    def write(self, output_mdp=None, skipempty=False):
        """
        Writes the MDP instance (the ordered dictionary) to an output MDP file.

        Parameters
        ----------
        output_mdp : str, Optional
            The file path of the output MDP file. The default is the filename the MDP instance was built from.
            If that if :code:`output_mdp` is not specified, the input MDP file will be overwritten.
        skipempty : bool, Optional
            Whether to skip empty values when writing the MDP file. If :code:`True`, any parameter lines from
            the output that contain empty values will be removed. The default is :code:`False`.
        """
        # The line 'if skipempty and (v == "" or v is None):' below could possibly incur FutureWarning
        warnings.simplefilter(action='ignore', category=FutureWarning)

        if output_mdp is None:
            output_mdp = self.input_mdp

        with open(output_mdp, "w") as mdp:
            for k, v in self.items():
                if k[0] == "B":  # blank line
                    mdp.write("\n")
                elif k[0] == "C":  # comment
                    mdp.write(f"; {v!s}\n")
                else:  # parameter = value
                    if skipempty and (v == "" or v is None):
                        continue
                    if isinstance(v, six.string_types) or not hasattr(v, "__iter__"):
                        mdp.write(f"{k!s} = {v!s}\n")
                    else:
                        mdp.write(f"{k} = {' '.join(map(str, v))}\n")
