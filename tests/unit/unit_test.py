'''
Created on 30 mars 2020

@author: cluzetb
'''
# content of test_sample.py
import os
from utilcrocO import unpack_conf
import pytest


def inc(x):
    return x + 1


def test_OMP_num_threads():
    """
    with a proper install, OMP_NUM_THREADS MUST be set to 1.
    """
    assert os.environ['OMP_NUM_THREADS'] == '1'


def test_unpack_conf():
    assert unpack_conf('aa,bb') == ['aa', 'bb']
    assert unpack_conf("'toto'") == "toto"
    assert unpack_conf("2") == 2
    assert unpack_conf("2.3") == 2.3
    assert unpack_conf('2.3,2.4') == [2.3, 2.4]
    assert unpack_conf("'aa','bb'") == ['aa', 'bb']
    assert unpack_conf("2016080106,2017080106") == ['2016080106', '2017080106']
    with pytest.raises(ValueError, match='"toto: badly formatted string.'):
        unpack_conf('"toto')
    with pytest.raises(TypeError, match="unconsistent types in list:*"):
        unpack_conf("aa,2013")
    assert unpack_conf("rangex(start:1 end:17)") == list(range(1, 18))
    assert unpack_conf('toto.nam') == 'toto.nam'
    assert unpack_conf('None') == None
