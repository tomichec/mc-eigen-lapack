/* Transition rates from Clancy, Rudy (2002) */
_VFUN(alpha11,3.802/(0.1027*exp(-V/17.0)+0.20*exp(-V/150.)))	/* C3->C2 (R->Q); IC3->IC2 (S->T) */
_VFUN(alpha12,3.802/(0.1027*exp(-V/15.0)+0.23*exp(-V/150.)))	/* C2->C1 (Q->P); IC2->IF (T->U) */
_VFUN(alpha13,3.802/(0.1027*exp(-V/12.0)+0.25*exp(-V/150.)))	/* C1->O (P->O) */
_VFUN(beta11,0.1917*exp(-V/20.3))				/* C2->C3 (Q->R); C2->C3 (T->S) */
_VFUN(beta12,0.20*exp(-(V-5)/20.3))				/* C1->C2 (P->Q); C1->C2 (U->T) */
_VFUN(beta13,0.22*exp(-(V-10)/20.3))				/* O->C1 (O->P) */
_VFUN(alpha3,3.7933e-7*exp(-V/7.7))				/* IF->C1 (U->P); IC2->C2 (T->Q); IC3->C3 (S->R) */
_VFUN(beta3,8.4e-3+2e-5*V)					/* C1->IF (P->U); C2->IC2 (Q->T); C3->IC3 (R->S) */
_VFUN(alpha2,9.178*exp(V/29.68))				/* O->IF (O->U) */
_VFUN(beta2,(alpha13*alpha2*alpha3)/(beta13*beta3))		/* IF->O (U->O) */
_VFUN(alpha4,alpha2/100.)					/* IF->IM1 (U->V) */
_VFUN(beta4,alpha3)						/* IM1->IF (V->U) */
_VFUN(alpha5,alpha2/(9.5e4))					/* IM1->IM2 (V->W) */
_VFUN(beta5,alpha3/50.)						/* IM2->IM1 (W->V) */

_RATE(O,C1,beta13,alpha13)
_RATE(C1,C2,beta12,alpha12)
_RATE(C2,C3,beta11,alpha11)
_RATE(C3,IC3,beta3,alpha3)
_RATE(IC3,IC2,alpha11,beta11)
_RATE(IC2,IC1,alpha12,beta12)
_RATE(IC1,IM1,alpha4,beta4)
_RATE(IM1,IM2,alpha5,beta5)
_RATE(IC2,C2,alpha3,beta3)
_RATE(IC1,C1,alpha3,beta3)
_RATE(IC1,O,beta2,alpha2)
