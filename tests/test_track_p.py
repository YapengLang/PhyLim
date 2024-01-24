# this test works for: 1)estimated tau by ens by bisect 2)grabing a random q, then renderring a p by expm. use allclose to check if this p
# is same as estimated by bisect.
from ylib.delta_col import tau_TX
from ylib.track_p import teleport as tp
import numpy
import pytest
from scipy.linalg import expm

Q = numpy.array(
    [
        [
            -0.09227284468529699,
            0.028045147578579357,
            0.032750077361831595,
            0.03147761974488603,
        ],
        [
            0.07371518704942183,
            -0.13795342549993928,
            0.015420545003243148,
            0.048817693447274296,
        ],
        [
            0.009679756935104084,
            0.009658426422327308,
            -0.052349144887941214,
            0.033010961530509815,
        ],
        [
            0.021873034165443374,
            0.028441815418702904,
            0.06130941737583104,
            -0.11162426695997732,
        ],
    ]
)

Q2 = numpy.array(
    [
        [
            -3.1913886545423408e-05,
            3.1913941083416741e-11,
            3.1913965917541626e-11,
            3.1913822717516412e-05,
        ],
        [
            6.3827645433246084e-03,
            -6.3827646071528818e-03,
            3.1914119575533403e-11,
            3.1914152850328605e-11,
        ],
        [
            3.1913965917541626e-11,
            3.1913822717516412e-05,
            -3.1913886545423408e-05,
            3.1913941083416741e-11,
        ],
        [
            3.1914119575533403e-11,
            3.1914152850328605e-11,
            6.3827645433246084e-03,
            -6.3827646071528818e-03,
        ],
    ]
)

P0 = numpy.array(
    [0.13050301681669385, 0.30086927600458385, 0.20384854752811415, 0.3647791596506081]
)

P02 = numpy.array(
    [0.1275121123107646, 0.3025142511046492, 0.18989463702953507, 0.3800789995550512]
)

taus_to_try = numpy.arange(1e-6, 50, 0.01)
ens_to_try = numpy.array(
    [
        0.5287255740798872,
        0.7080280055427696,
        0.7208769669532618,
        0.7590550112580364,
        0.8234346574159457,
        0.829388321451194,
        0.861349466906161,
        0.8735144138940558,
        0.8811377961243829,
        0.8855647004227656,
        0.8875941839556121,
        0.8925237336662164,
        0.8978888293078392,
        0.9097734392674511,
        0.9155354500809586,
        0.9211312026787605,
        0.9259777424055387,
        0.9333316561225572,
        0.9500959156210316,
        0.9516972682232615,
        0.9593547067891117,
        0.9652796774488029,
        0.9850819143667263,
        1.0400449073096392,
        1.0511716791420471,
    ]
)


@pytest.mark.parametrize("tau", taus_to_try)
def test_tau_tx(tau):
    where = tau_TX(p0=P0, q=Q, tau=tau)
    assert numpy.allclose(tau, tau_TX(p0=P0, q=Q, ens=where))

    where = tau_TX(p0=P02, q=Q2, tau=tau)
    assert numpy.allclose(tau, tau_TX(p0=P02, q=Q2, ens=where))


@pytest.mark.parametrize("tau", taus_to_try)
def test_est_p(tau):
    expected = expm(Q * tau)

    where_ens = tau_TX(p0=P0, q=Q, tau=tau)
    est_tau = tau_TX(p0=P0, q=Q, ens=where_ens)

    assert numpy.allclose(expected, expm(Q * est_tau))

    expected = expm(Q2 * tau)

    where_ens = tau_TX(p0=P02, q=Q2, tau=tau)
    est_tau = tau_TX(p0=P02, q=Q2, ens=where_ens)

    assert numpy.allclose(expected, expm(Q2 * est_tau))


@pytest.mark.parametrize("ens", ens_to_try)
def test_track_p(ens):
    tp(M="ssGN", want_chainsaw=True, where_ens=ens)
