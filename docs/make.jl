

#push!(LOAD_PATH,"../src/")

using Documenter, ANISPROU


makedocs(sitename="ANISPROU",
         modules = [ANISPROU],
         authors = "Andrea Zunino",
         format = Documenter.HTML(prettyurls=get(ENV,"CI",nothing)=="true")  )


deploydocs(
    repo="github.com/inverseproblem/ANISPROU.jl.git",
    devbranch = "main",
    branch = "gh-pages" 
)
