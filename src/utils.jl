

####################################################

"""
$(TYPEDSIGNATURES)

Least squares linear regression.
"""
function lssqregr(points::Array{<:Real,2})
    
    ## G = [protc ones(length(protc))]
    ## dobs = sdsc
    G = [points[:,1] ones(size(points,1))]
    dobs = points[:,2]

    ## overdetermined least squares
    ## mpost = (G'*G)*G'*dobs
    ## (G'*G) * h = G'
    lstmp = (G'*G) \ G'
    mpost = lstmp * dobs

    angcoe = mpost[1]
    intercept = mpost[2]

    return angcoe,intercept
end

####################################################
