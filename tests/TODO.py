import pytest


# HIGH LEVEL FUNCTIONS:
#   found in xtea.locate_insertions
#      These functions write out files, check output of files for matching (off by one okay)

# get_clip_sites
# get_disc_sites
# filter_csn


####### samples
def f():
    raise SystemExit(1)


def test_mytest():
    with pytest.raises(SystemExit):
        f()