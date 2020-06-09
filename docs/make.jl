

push!(LOAD_PATH,"../src/")

using Documenter, ANISPROU


makedocs(modules = [ANISPROU],
         sitename="ANISPROU",
         format = Documenter.HTML(prettyurls=get(ENV,"CI",nothing)=="true")  )

