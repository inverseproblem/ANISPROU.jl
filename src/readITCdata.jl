

###########################################################

"""
$(TYPEDSIGNATURES)

Real from ASCII files all observed/measured data from a set of experiments.
It can read data for one or more proteins.
The entalphy values are *scaled* by a factor given by the argument `scalfactor`, 
 which defaults to 0.004184 (Cal/mol to kJ/mol).

# Arguments 

-`inpdir`: directory containing the input data
-`proteinnames`: array of strings containing the names of proteins
-`scalfactor`=0.004184: scaling factor for enthalpy, defaults 
                        to 0.004184 (Cal/mol to kJ/mol)
-`discninitrows`=0: number of initial rows of the data set to discard. This is used
                    to remove some initial data often affected by strong 
                    instrument noise which could bias the fitting process.

"""
function readallexperiments(inpdir::String,proteinnames::Vector{String} ;
                            scalfactor::Float64=0.004184 ,
                            discninitrows::Integer=0)
    ## swith units, from Cal/mol to kJ/mol by multiplying by 0.004184
    @assert discninitrows>=0
    starow = discninitrows+1
    #@show discninitrows

    println("Reading data from directory: $inpdir")
    println(" Scaling enthalpy values by a factor $scalfactor")

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
        # @show fllist
        # @show occursin(proteinnames[p],fllist[1])
        # @show endswith(fllist[1],".DAT")
        # select only files starting with
        checkfile(x) = (occursin(proteinnames[p],x) && endswith(x,".DAT"))
        fllist = filter(checkfile,fllist)
        println(" file list: \n$(fllist)")

        startind = 1
        for fl in fllist

            seflname = joinpath(inpdir,proteinnames[p]*"/"*fl)
            csds1,cpro1,ceout1 = readsingleexperiment(seflname)

            # discard some initial rows (strongly affected by instrument error)
            csds1 = csds1[starow:end] 
            cpro1 = cpro1[starow:end]
            ceout1 = ceout1[starow:end]

            # add to dictionary
            lendata = length(ceout1)
            push!(sdscon,csds1...)
            push!(procon,cpro1...)
            push!(enout,ceout1...)
            # indices of individual experiments
            push!(idxdata,startind:(startind+lendata-1))
            startind += lendata
        end

        ## scale enthalpy
        enout .*= scalfactor

    end

    return data
end

###########################################################

"""
$(TYPEDSIGNATURES)

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

    # REMARK:
    # For DH, INJV, XMt and NDH an initial row with a zero should be added.
    # Xt and Mt are initially one row longer than the other columns, and
    #   this is corrected for by the additional zeros.

    data,header = readdlm(singlefl,header=true)
    header = vec(header)

    idxXt = findall(x->x=="Xt",header)[1]
    idxMt = findall(x->x=="Mt",header)[1]
    @assert idxXt<idxMt

    ## DH	INJV	Xt	Mt	XMt	NDH
    # the last row contains only Xt and Mt at the beginning...
    data2 = data[1:end-1,:]
    lasttXtMt = data[end,:]
   
    idx=findall(x->x=="Xt",header)[1]
    # Xt in the last row of data is in the second column
    sdscon = convert.(Float64,[data2[:,idx]; lasttXtMt[2]])

    idx=findall(x->x=="Mt",header)[1]
    # Xt in the last row of data is in the third column
    procon = convert.(Float64,[data2[:,idx]; lasttXtMt[3]])

    idx=findall(x->x=="NDH",header)[1]
    enout = convert.(Float64,[0.0; data2[:,idx]])

    return sdscon,procon,enout
end

###########################################################
