r"""
Copyright (c) 2023 J. K. Eriksen, L. Panny, J. Sotáková, M. Veroni

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import collections

from sage.all import *
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
#from sage.misc.verbose import verbose


from .klpt import KLPT_Context, DecompAlphaN

#####################################
#                                   #
#            Main/Wrapper           #
#                                   #
#####################################

def constructive_deuring(I, E0, iota, variant=None):
    r"""
    Given a curve E0, where End(E0) is isomorphic to O0 (a special
    p-extremal maximal order in the quaternion algebra B = <1,i,j,k>),
    the endomorphism iota corresponding to i, and a left O0 ideal I,
    outputs a curve E_I, whose endomorphism ring is isomorphic to the right order of I,
    and an isogeny phi_J : E_0 -> E_I of smooth norm.
    """
    O0 = I.left_order()
    KLPT_ctx = KLPT_Context(O0.quaternion_algebra())
    if not sage.misc.banner.require_version(9,8):
        EllipticCurveHom_composite.make_default()  #remove once everyone runs Sage >= 9.8

    J, facToExt, T, S, f = KLPT_ctx.KLPT(I, variant=variant)
    #print(f'norm(J) = {factor(J.norm())}')
    #verbose(f'{J=}')
    assert gcd(J.left_order().unit_ideal().basis_matrix().solve_left(J.basis_matrix()).list()) == 1, 'non-cyclic ideal'

    grouped = collections.defaultdict(list)
    for le,k in sorted(facToExt.items()):
        grouped[k].append(le)
    #print('Torsion by extension degree:')
    #for k,les in sorted(grouped.items()):
        #print(f'  {k:4}:  ', end='')
        #for j,le in enumerate(les):
            #print((', ' if j else '') + str(factor(le)), end='')
        #print()

    Deuring_ctx = Deuring_Context(O0, E0, iota, facToExt, T, S, f, J)

    if S:
        #print("Using fancier translation")
        #print(f"Translating {f} ideals of norm dividing S = {factor(S)}, by coprime T = {factor(T)}")
        phi_I = Deuring_ctx.IdealToIsogeny_S(J)

    else:
        #print("Using one-shot ideal-isogeny translation")
        phi_I = Deuring_ctx.IdealToIsogeny(J)

    return phi_I.codomain(), phi_I, Deuring_ctx


#####################################
#                                   #
#     Ideal-Isogeny translation     #
#                                   #
#####################################

class TwistedCurveWithEndomorphism:
    def __init__(self, E, D, iota_maps):
        _,_,_,a,b = E.a_invariants()
        if E != EllipticCurve([a,b]):
            raise NotImplementedError

        self.E = EllipticCurve([a/D**2, b/D**3])
        self.D = D

        (fnum,fden), (gnum,gden) = iota_maps
        x,y = Sequence([fnum,fden,gnum,gden]).universe().gens()
        self.f = fnum(x=D*x), fden(x=D*x)*D
        assert gnum.coefficient({y:1})*y == gnum
        assert set(gden.variables()) <= {x}
        self.g = gnum.coefficient({y:1})(x=D*x)*y, gden(x=D*x)

        self.p = E.base_field().characteristic()
        self.u = D**(self.p//2)

    def eval_i(self, pt):
        assert pt in self.E
        if not pt:
            return pt
        x,y = pt.xy()
        return self.E(self.f[0](x,y)/self.f[1](x,y), self.g[0](x,y)/self.g[1](x,y))

    def eval_j(self, pt):
        assert pt in self.E
        if not pt:
            return pt
        x,y = pt.xy()
        return self.E(x**self.p * self.u**2, y**self.p * self.u**3)

class Deuring_Context:
    r"""
    Helper to setup parameters for computing the Deuring correspondence.
    """
    def __init__(self, O0, E0, iota, facToExt, T, S=None, f=None, J=None):
        self.O0 = O0
        self.E0 = E0
        F2 = E0.base_field()
        assert F2.degree() == 2
        self.p = F2.characteristic()
        self.iota = iota
        self.facToExt = facToExt
        self.T = T
        self.S = S
        self.f = f
        self.J = J

    @cached_method
    def extE0(self, extdeg, twist):
        r"""
        Returns either a base-change of E_0 over F_{p^2k},
        or a quadratic twist \widetilde{E_0} of E_0 over F_{p^2k}.
        The twist has its own isomorphism \widetilde{E_0} -> O_0
        """
        Fbig,emb = self.E0.base_field().extension(extdeg,'A', map=True)
        E0big = self.E0.change_ring(Fbig)
        if twist:
            while True:
                D = Fbig.random_element()
                if not D.is_square():
                    break
        else:
            D = Fbig.one()

        iota_maps = tuple(tuple(g.map_coefficients(emb)
                                for g in (f.numerator(), f.denominator()))
                          for f in self.iota.rational_maps())

        return TwistedCurveWithEndomorphism(E0big, D, iota_maps), emb

    def evalIdealElt(self, EE, a, Q):
        r"""
        Given an endomorphism a and a point Q of
        order coprime to the denominator of a,
        returns a(Q)
        """
        d = lcm(c.denominator() for c in a)
        iQ = EE.eval_i(Q)
        jQ = EE.eval_j(Q)
        kQ = EE.eval_i(jQ)
        coeffs = [coeff % Q.order() for coeff in a]
        aQ = coeffs[0]*Q + coeffs[1]*iQ + coeffs[2]*jQ + coeffs[3]*kQ
        return aQ

    def _precompute_torsion_bases(self, Eles):
        from os import cpu_count
        for (Ele,_), (EE,emb, P,Q) in parallel(cpu_count())(self.genTorsionBasis)(Eles):
            self.genTorsionBasis.set_cache((EE,emb, P,Q), *Ele)
            assert emb.domain() is self.E0.base_ring()
            assert emb.codomain() is P.base_ring() is Q.base_ring()
            if (emb0 := emb.codomain().coerce_map_from(emb.domain())) is None:
                emb.codomain().register_coercion(emb)
            assert emb.codomain().coerce_map_from(emb.domain()) == emb

    def IdealToIsogenyGens(self, I, specificTorsion=0):
        r"""
        Given a left O0-ideal I, returns {P_1, .., P_n} such that ker phi_I = <P_1, .., P_n>
        """
        kerGens = []
        a = DecompAlphaN(I)
        d = lcm(c.denominator() for c in a)
        N = ZZ(I.norm()).gcd(specificTorsion) if specificTorsion else ZZ(I.norm())

        self._precompute_torsion_bases([(self.E0, l, e+d.valuation(l)) for l,e in N.factor()])

        for (l,e) in N.factor():
            lval = d.valuation(l)  # point divisions
            EE,emb, P,Q = self.genTorsionBasis(self.E0, l, e+lval)
            assert P in EE.E and Q in EE.E
            assert EE.E.base_ring().has_coerce_map_from(self.E0.base_ring())
            R = self.evalIdealElt(EE, l**lval * a.conjugate(), P)
            if not l**(e-1) * R:
                R = self.evalIdealElt(EE, l**lval * a.conjugate(), Q)
            kerGens.append((xPoint(R.xy()[0]*EE.D, self.E0.change_ring(emb)), (l,e)))
        return kerGens

    def IdealToIsogeny(self, I, specificTorsion=0):
        r"""
        Given a left O0-ideal I, returns the corresponding isogeny phi_I
        """
        kerGens = self.IdealToIsogenyGens(I, specificTorsion=specificTorsion)
        #verbose("Evaluating isogeny...")
        return chain_iso(kerGens, self.E0)


    #####################################
    #                                   #
    # Sliding translation (Appendix A)  #
    #                                   #
    #####################################

    def dualIsoGens(self, phi, degree):
        r"""
        Given an isogeny phi and its degree,
        outputs a generating set of the kernel of the dual isogeny
        """
        facDeg = factor(degree)
        kerGens = []
        for (l,e) in facDeg:
            P, Q = self.genTorsionBasis(phi.domain(), l, e, xOnly=True)
            #verbose(f"Finding kernel generator of order {l}^{e} over F_p^{P.curve.base_field().degree()}")
            R = P.push(phi)
            if not R.mul(l**(e-1)):
                R = Q.push(phi)
            kerGens.append((R, (l,e)))
        return kerGens

    def dualIso(self, phi, degree):
        r"""
        Given an isogeny phi and its degree,
        returns the dual of phi up to automorphisms
        """
        kerGens = self.dualIsoGens(phi, degree)
        phihat = chain_iso(kerGens, phi.codomain())
        varphi = phihat.codomain().isomorphism_to(phi.domain()) 
        return varphi * phihat

    def IdealToIsogenyCoprime(self, J, K, phi_K, alpha):
        r"""
        Given left O0-ideals J, K of coprime norms, phi_K and alpha such that
        J = K*\bar{alpha}/nrd(K), outputs phi_J
        """
        H1 = J + self.O0*self.T
        H2 = self.O0*alpha + self.O0*(J.norm()/H1.norm())

        #verbose("--> Translating H1...")
        phi_H1 = self.IdealToIsogeny(H1)

        #verbose("--> Translating H2...")
        ker_psi = self.IdealToIsogenyGens(H2)
        ker_psi = [(R.push(phi_K), deg) for (R, deg) in ker_psi]
        psi = chain_iso(ker_psi, phi_K.codomain())
        assert psi.codomain().j_invariant() == phi_H1.codomain().j_invariant()

        #verbose("--> Moving precomputed basis")
        precompbasis = []
        for l, e in factor(H2.norm()):
            assert not l.divides(phi_K.degree())
            basis = self.genTorsionBasis(phi_K.domain(), l, e, xOnly=True)
            impts = tuple(T.push(phi_K) for T in basis)
            self.genTorsionBasis.set_cache(impts, phi_K.codomain(), l, e, xOnly=True)
        #verbose("--> Computing the dual")
        psihat = self.dualIso(psi, psi.degree())
        varphi = phi_H1.codomain().isomorphism_to(psihat.domain())
        assert psihat.codomain() == phi_K.codomain()
        return psihat * varphi * phi_H1

    def IdealSlide(self, I, J, phi_J):
        r"""
        Given left O0-ideals I, J, phi_J, where I and J have coprime norm,
        such that I = J*Ii for some ideal Ii of norm S, outputs phi_Ii
        """
        assert self.O0 == I.left_order()
        I1 = I + self.O0*self.S

        #print(f"Translating Ii of norm {factor(I1.norm())} to kernel")
        ker_I1 = self.IdealToIsogenyGens(I1)
        #print("Pushing ker phi_I1 through phi_J")

        ker_phi1 = [(R.push(phi_J), deg) for (R, deg) in ker_I1]
        return chain_iso(ker_phi1, phi_J.codomain())

    def IdealFiltration(self, I):
        r"""
        Given an ideal I of norm dividing S^f, returns f ideals I_i,
        each of norm dividing S, such that I = I_1*I_2*...*I_f
        """
        ideals = []
        assert self.O0 == I.left_order()
        Ibuild = self.O0.unit_ideal()
        for i in range(1, self.f+1):
            I_i = Ibuild.conjugate()*I*(1/Ibuild.norm()) + Ibuild.right_order()*self.S
            Ibuild = Ibuild * I_i
            #verbose(f"I_{i-1} has norm {factor(I_i.norm())}")
            ideals.append(I_i)
        assert Ibuild == I
        return ideals

    @cached_method
    def genTorsionBasis(self, E, l, e, xOnly=False):
        ### is there anything to be gained by using a more clever basis
        ### such that we do know the action of the endo ring on it?
        ### eg SIDH choice P, iota(P)
        r"""
        Given a curve E, prime l, exponent e, output generators (P, Q) of E[l^e]
        """
        t = l**e
        extdeg = self.facToExt[t]
        twist = not t.divides(self.p**extdeg - (-1)**extdeg)
        #verbose(f"Generating torsion basis for E[{l}^{e}] over F_p^{2*extdeg}" + (" on quadratic twist" if twist else ""))

        F = E.base_field()
        EE,emb = self.extE0(extdeg, twist)

        order = self.p**extdeg - (-1 if twist else +1) * (-1)**extdeg
        EE.E.set_order(order**2, num_checks=0)
        assert t.divides(order)

        cof = order.prime_to_m_part(l)

        def rpt():
            while True:
                T = cof * EE.E.random_point()
                Tl = l**(e-1)*T
                if Tl: break
            U = l*Tl
            while U:
                Tl = U
                U *= l
                T *= l
#            assert l**(e-1)*T and not l**e*T
#            assert Tl and not l*Tl
            T.set_order(l**e)
            Tl.set_order(l)
            return T, Tl

        P,Pl = rpt()
        Q,Ql = rpt()
        while Pl.weil_pairing(Ql,l).is_one():
            Q,Ql = rpt()

        if xOnly:
            P = xPoint(P.xy()[0]*EE.D, self.E0.change_ring(emb))
            Q = xPoint(Q.xy()[0]*EE.D, self.E0.change_ring(emb))
            return P, Q

        return EE,emb, P,Q

    def IdealToIsogeny_S(self, I):
        r"""
        Given a left O0-ideal I of norm S^f, outputs phi_I
        """
        ctx = KLPT_Context(I.quaternion_algebra())

        #verbose("----> Creating Ideal filtration")
        Ifilt = self.IdealFiltration(I)
        K = Ifilt[0]

        #verbose(f"----> Translating K (of norm {factor(K.norm())}):")
        phi_K = self.IdealToIsogeny(K)

        #verbose("----> Finding J equivalent (and of coprime norm) to K")
        beta, J, _, _, _, _ = ctx.KLPT(K, T=self.T**2, returnElem=True)

        #verbose("----> Finding phi_J through IdealToIsogenyCoprime")
        phi_J = self.IdealToIsogenyCoprime(J, K, phi_K, beta)

        phiouts = []
        for index in range(1, len(Ifilt)):
            #print(f"----> Passing part {index} out of {len(Ifilt) - 1}")
            Ii = multiplyIdeals(J, Ifilt[index], beta=beta.conjugate())
            phi_I1 = self.IdealSlide(Ii, J, phi_J)
            phiouts.append(phi_I1)

            if index == len(Ifilt) - 1:
                break

            K = K*Ifilt[index]
            #print(f'norm(K) = {factor(K.norm())}')
            phi_K = phi_I1 * phi_K

            #verbose("----> Finding J equivalent (and of coprime norm) to K")
            beta, J, _, _, _, _ = ctx.KLPT(K, T=self.T**2, returnElem=True)

            #verbose("----> Finding phi_J through IdealToIsogenyCoprime")
            phi_J = self.IdealToIsogenyCoprime(J, K, phi_K, beta)

        return EllipticCurveHom_composite.from_factors(phiouts)


#####################################
#                                   #
#       Ideal helper functions      #
#                                   #
#####################################

def multiplyIdeals(I, J, beta=None):
    r"""
    Computes the product I*J of two quaternion ideals I, J.
    If the right order of I and the left order of J are of
    same type, but not equal, an isomorphism orders must be given.

    By Skolem-Noether, this isomorphism is given by conjugation by
    a quaternion beta.
    """
    if I.right_order() != J.left_order():
        assert beta, "Must provide automorphism"
        J = ~beta*J*beta
        assert I.right_order() == J.left_order(), "Orders does still not match"
    return I*J


#####################################
#                                   #
#           EC Functions            #
#                                   #
#####################################

from .xonly import xPoint, xISOG

def chain_iso(kernelPointsIn, E):
    r"""
    Given points {P_1, ..., P_n} on curve E, outputs the isogeny
    whose kernel is generated by {P_1, ..., P_n}
    """
    Ei = E
    kernelPoints = kernelPointsIn[:]

    philist = []
    while kernelPoints:
        kernelPoints.sort(key = lambda tup: -tup[1][0])  # start with small degrees
        Ri, (l,e) = kernelPoints.pop()
        #verbose(f'Computing {l}-isogeny step, pushing {len(kernelPoints)+(e>1)} points')
        Ki = Ri.mul(l**(e-1))
        phi = xISOG(Ei, Ki.X, l)
        Ei = phi.codomain()

        if e > 1:
            kernelPoints.append((Ri, (l, e-1)))

        from os import cpu_count
        chunks = [[] for _ in range(cpu_count())]
        for i,tup in enumerate(kernelPoints):
            chunks[i % len(chunks)].append(tup)
        def push(tups):
            return [(P.push(phi), order) for P, order in tups]
        if len(kernelPoints) > 1:
            kernelPoints = [tup for _,tups in parallel(cpu_count())(push)(chunks) for tup in tups]
        else:
            kernelPoints = push(kernelPoints)

        philist.append(phi)

    return EllipticCurveHom_composite.from_factors(philist)

