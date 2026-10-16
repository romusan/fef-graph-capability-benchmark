"""Checks against exact distributions, identities and published archived counts."""
import numpy as np
from scipy import stats
import reanalisis_estadistico as a

def test_signed_rank_extremes_and_ties():
    assert a.rank_biserial([1,2,3])==1
    assert a.rank_biserial([-1,-2,-3])==-1
    assert a.rank_biserial([1,-1,0])==0
    assert a.rank_biserial([0,0])==0

def test_exact_mcnemar_and_orientation():
    p,odds,ci=a.mcnemar_exact(0,6,0,24)
    assert p==.03125 and odds==0
    p,odds,ci=a.mcnemar_exact(14,11,5,0)
    assert abs(p-.210113525390625)<1e-14
    assert odds==5/11

def test_exact_power_reduces_to_binomial_tail():
    # All discordance points in one direction. Rejection iff D >= 6.
    expected=stats.binom.sf(5,30,.2)
    assert abs(a.power_mcnemar(30,.2,.2)-expected)<1e-14
    assert a.emd_mcnemar(30,.2) is None
    assert a.power_mcnemar(30,.5,0)<=.05

def test_holm_monotone_and_qd_schema():
    np.testing.assert_allclose(a.holm([.01,.02,.8]),[.03,.04,.8])
    assert a.resultados['qd_arms']['Estrict']['strict_seeds']==9
    assert a.resultados['qd_arms']['Bmatched']['strict_seeds']==10
    assert a.resultados['qd_arms']['Estrict_vs_Bmatched_pareado']['n_comunes']==15

def test_paired_effect_and_holm_archived_regression():
    p=a.resultados['packaging']
    assert abs(p['wattI']['rms']['rank_biserial']-.5010752688172043)<1e-12
    assert abs(p['gearedFiveBar']['rms']['rank_biserial']+.2)<1e-12
    assert abs(p['holm_packaging']['wattI:rms']['p_holm']-.09276001900434494)<1e-12
