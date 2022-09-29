#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def attachment_rate_cte(rho,km,k0):
    kp = (km*rho*(1-rho)-k0*rho**2)/(1-2*rho)**2
    return kp




