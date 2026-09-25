"""ddsim.tools has no functional test coverage (see docs/architecture.md) -- at
least confirm its modules import cleanly. This is the only thing that would
have caught the ddsim.tools.PDF/twins relative-import bug (`from . import X`
instead of `from .. import X`, since tools/ is one level deeper than the
modules it imports) introduced when the package was restructured into
src/ddsim/.
"""


def test_tools_modules_import():
    from ddsim.tools import PDF, twins  # noqa: F401
