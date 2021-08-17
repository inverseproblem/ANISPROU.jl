

####################################################

"""
$(TYPEDSIGNATURES)

Least squares linear regression in its simplest form.
"""
function lssqregr(points::Array{<:Real,2})
    
    ## G = [protc ones(length(protc))]
    ## dobs = sdsc
    G = [points[:,2] ones(size(points,1))]
    dobs = points[:,1]

    ## overdetermined least squares
    ## mpost = (G'*G)*G'*dobs
    ## (G'*G) * h = G'
    lstmp = (G'*G) \ G'
    mpost = lstmp * dobs

    angcoe = mpost[1]
    intercept = mpost[2]

    # residuals
    r = (G * mpost) .- dobs
    sse = sum(r.^2)
    if length(r)>2
        resstderr = sqrt(sse/(length(dobs)-2))
    else
        resstderr = sqrt(sse/(length(dobs)))
    end
    
    return angcoe,intercept,resstderr
end

####################################################

# """
# $(TYPEDSIGNATURES)

# Least squares linear regression including a covariance 
#  matrix on the observed/measured data.

# Cd is the covariance matrix on the observed data.
# """
# function lssqregr(points::Array{<:Real,2},Cd::AbstractMatrix{<:Real})

#     # take the inverse of the covariance matrix
#     invCd = inv( Cd )

#     ## G = [protc ones(length(protc))]
#     ## dobs = sdsc
#     G = [points[:,2] ones(size(points,1))]
#     dobs = points[:,1]
    
#     ## overdetermined least squares
#     ## mpost = (G'*G)*G'*dobs
#     ## (G'*G) * h = G'
#     lstmp = (G'*invCd*G) \ (G'*invCd)
#     mpost = lstmp * dobs

#     # see Tarantola 2005, page 67, eq. 3.40-3.44
#     # posterior covariance on model parameters
#     Cmpost = inv( transpose(G) * invCd * G )
#     # # posterior covariance on observed data
#     # Cdpost = G * Cmpost * transpose(G)

#     angcoe = mpost[1]
#     intercept = mpost[2]

#     # posterior variance
#     #varangcoe,varintercept = Cmpost[1,1],Cmpost[2,2]

#     # residuals
#     r = (G * mpost) .- dobs
#     sse = sum(r.^2)
#     if length(r)>2
#         resstderr = sqrt(sse/(length(dobs)-2))
#     else
#         resstderr = sqrt(sse/(length(dobs)))
#     end
    

#     @show G
#     @show dobs
    
#     return angcoe,intercept,resstderr #varangcoe,varintercept
# end

# ####################################################
