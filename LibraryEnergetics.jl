using LinearAlgebra: norm

@inline perturbation_norm(field, bkgd_field) = norm(field - bkgd_field)