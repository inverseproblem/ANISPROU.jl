

###########################################################

"""
$TYPEDSIGNATURE)

Real from ASCII files all observed/measured data from a set of experiments.
It can read data for one or more proteins.
"""
function readallexperiments(inpdir::String,proteinnames::Vector{String})

    numprot = length(proteinnames)

    data = Dict()
    

    for p=1:numprot

        data[proteinnames[p]] = Dict()
        data[proteinnames[p]]["sdscon"] = Float64[]
        data[proteinnames[p]]["procon"] = Float64[]
        data[proteinnames[p]]["enout"] = Float64[]
        data[proteinnames[p]]["idxdata"] = UnitRange{Int64}[]

        sdscon = data[proteinnames[p]]["sdscon"]
        procon = data[proteinnames[p]]["procon"]
        enout = data[proteinnames[p]]["enout"]
        idxdata = data[proteinnames[p]]["idxdata"] 

        # get a list of files in inpdir
        fllist = readdir(joinpath(inpdir,proteinnames[p]))
        # select only files starting with
        checkfile(x) = (startswith(x,proteinnames[p]) && endswith(x,"uM.txt"))
        fllist = filter(checkfile,fllist)

        startind = 1
        for fl in fllist
            csds1,cpro1,ceout1 = readsingleexperiment(joinpath(inpdir,proteinnames[p]*"/"*fl))
            lendata = length(ceout1)
            push!(sdscon,csds1...)
            push!(procon,cpro1...)
            push!(enout,ceout1...)
            @show startind:(startind+lendata-1),lendata
            push!(idxdata,startind:(startind+lendata-1))
            startind += lendata
        end

    end

    return data
end

###########################################################

"""
$TYPEDSIGNATURE)

Real from ASCII file the observed/measured data from an experiment.
"""
function readsingleexperiment(singlefl::String)
    ## Only column 3 (Xt), 4 (Mt) and 6 (NDH) are relevant
    #
    # Xt is the concentration of SDS in the cell (in mM). This
    # is initially 0, and increases during the experiment.
    #
    # Mt is the protein concentration. This is initially the
    # value stated in the header of the file (in mM), and
    # decreases slightly during the experiment (as it is
    # diluted with the SDS).
    #
    # NDH is the energy output of that single injection in J/mol..
    # 
    
    data,header = readdlm(singlefl,header=true)
    header = vec(header)

    idx=findall(x->x=="Xt",header)[1]
    sdscon = data[:,idx]

    idx=findall(x->x=="Mt",header)[1]
    procont = data[:,idx]

    idx=findall(x->x=="NDH",header)[1]
    enout = data[:,idx]

    return sdscon,procont,enout
end

###########################################################
