module Conf


using ComfyCommons: info, warn, Git, Yaml
import Util, Data


"""
    all()

Create all job configurations.
"""
all() = (stepsize();  comparison())


"""
    stepsize([metaconfig = "conf/job/meta/stepsize.yml"])

Generate a set of job configurations from the given meta-configuration file.
"""
function stepsize(metaconfig::String="conf/job/meta/stepsize.yml")
    meta = load_file(metaconfig) # read config
    
    # check for recent changes
    githaschanges = Git.haschanges("src", metaconfig)
    if (githaschanges) warn("Uncommited changes may affect configurations") end
    
    # expand job configurations
    for meta_dataset in Yaml.expand(meta, "dataset"),
                 job in Yaml.expand(meta_dataset, ["dataset", "skconfig"])
        
        # interpolate
        interpolate!(job, "configfile", id = job["dataset"]["id"], skconfig = job["dataset"]["skconfig"])
        interpolate!(job, "outfile",    id = job["dataset"]["id"], skconfig = job["dataset"]["skconfig"])
        
        # expand each stepsize method and recombine
        job["experiments"] = vcat(map(
            exp -> begin
                exp_method = exp["stepsize"]["method"]
                if exp_method == "constant"
                    Yaml.expand(exp, ["stepsize", "alpha"])
                elseif in(exp_method, ["decay_mul", "decay_exp"])
                    Yaml.expand(exp, ["stepsize", "eta"])
                elseif exp_method == "run"
                    Yaml.expand(exp, ["stepsize", "tau"])
                else # e.g., optimal
                    exp
                end
            end, job["experiments"])...)
        
        # write job to file
        info("Writing configuration of $(length(job["experiments"])) experiments",
             " to $(job["configfile"])")
        prefix = """
        # 
        # $(job["configfile"])
        # 
        # Automatically generated with Conf.stepsize(\"$metaconfig\").
        # 
        # Git commit: $(Git.commithash())
        # Git origin: $(Git.remoteurl())
        # Git uncommited changes: $githaschanges
        # 
        """
        write_file(job["configfile"], job, prefix)
        
    end
end


"""
    comparison([metaconfig = "conf/job/meta/comparison.yml"])

Generate a set of job configurations from the given meta-configuration file.
"""
function comparison(metaconfig::String="conf/job/meta/comparison.yml")
    
    # read config and check for recent changes
    meta = load_file(metaconfig)
    githaschanges = Git.haschanges("src", metaconfig)
    if (githaschanges) warn("Uncommited changes may affect configurations") end
    
    # expand configuration
    for meta_dataset in Yaml.expand(meta, "dataset"),
                 job in Yaml.expand(meta_dataset, ["dataset", "f_train"])
            
        # interpolate
        dataset_id = job["dataset"]["id"]
        f_train    = job["dataset"]["f_train"]
        interpolate!(job, "configfile", id = dataset_id, ftrain = f_train)
        interpolate!(job, "outfile",    id = dataset_id, ftrain = f_train)
        
        # expand and filter experiments
        dim_f = length(Data.bins(Data.dataset(dataset_id, readdata = false))) # dimension of f
        job["experiments"] = filter(exp -> begin
            
            # smoothing order has to be smaller than dim_f
            if exp["method"] == "ibu" && exp["smoothing"]["method"] == "polynomial"
                exp["smoothing"]["order"] <=  dim_f - 1
            
            # num_clusters should be greater than dim_f
            elseif exp["method"] == "run" && exp["discretization"]["method"] == "expand"
                exp["discretization"]["factor"] * dim_f <= exp["num_clusters"]
            
            # omit redundant config of factor 1
            elseif exp["method"] == "run" && exp["discretization"]["method"] == "reduce"
                exp["discretization"]["factor"] != 1
            
            # keep the rest
            else  true  end
            
        end, vcat(map(exp -> begin
            
            # stepsize expansions for DSEA
            if exp["method"] == "dsea"
                exp_method = exp["stepsize"]["method"] # stepsize method
                if exp_method == "constant"
                    Yaml.expand(exp, ["stepsize", "alpha"])
                elseif in(exp_method, ["decay_mul", "decay_exp"])
                    Yaml.expand(exp, ["stepsize", "eta"])
                elseif exp_method == "run"
                    Yaml.expand(exp, ["stepsize", "tau"])
                else
                    throw(ArgumentError("Illegal stepsize method $(exp_method)"))
                end
            
            # expand order and num_clusters for IBU
            elseif exp["method"] == "ibu"
                exp_method = exp["smoothing"]["method"] # smoothing method
                if exp_method == "polynomial"
                    Yaml.expand(exp, ["smoothing", "order"], "num_clusters")
                elseif exp_method == "none"
                    Yaml.expand(exp, "num_clusters")
                else
                    throw(ArgumentError("Illegal smoothing method $(exp_method)"))
                end
            
            # expand method and factor of discretization and num_clusters for RUN
            elseif exp["method"] == "run"
                Yaml.expand(exp, ["discretization", "method"],
                                 ["discretization", "factor"],
                                 "num_clusters")
            
            # other methods not supported
            else  throw(ArgumentError("Illegal method $(exp["method"])"))  end
            
        end, job["experiments"])...))
        
        
        # informative name of strategy requires interpolation
        for exp in job["experiments"]
            name = exp["name"]
            
            if exp["method"] == "dsea"
                stepsize_conf = exp["stepsize"]
                exp_method = stepsize_conf["method"] # stepsize method
                if exp_method == "constant"
                    name = replace(name, "\$alpha", stepsize_conf["alpha"])
                elseif in(exp_method, ["decay_mul", "decay_exp"])
                    name = replace(name, "\$eta", stepsize_conf["eta"])
                elseif exp_method == "run"
                    name = replace(name, "\$tau", Util.latex_e(stepsize_conf["tau"], dollars = false))
                end
                
            elseif exp["method"] == "ibu"
                name = replace(name, "\$num_clusters", exp["num_clusters"])
                if exp["smoothing"]["method"] == "polynomial"
                    name = replace(name, "\$order", exp["smoothing"]["order"])
                end
            
            elseif exp["method"] == "run"
                name = replace(name, "\$num_clusters", exp["num_clusters"])
                if exp["discretization"]["factor"] == 1
                    name = replace(name, "\$method with factor \$factor", "no regularization")
                end
                name = replace(name, "\$method", exp["discretization"]["method"])
                name = replace(name, "\$factor", exp["discretization"]["factor"])
            
            end
            exp["name"] = name # replace with interpolation
        end
        
        
        # write job to file
        info("Writing configuration of $(length(job["experiments"])) experiments",
             " to $(job["configfile"])")
        prefix = """
        # 
        # $(job["configfile"])
        # 
        # Automatically generated with Conf.expandstepsize(\"$metaconfig\").
        # 
        # Git commit: $(Git.commithash())
        # Git origin: $(Git.remoteurl())
        # Git uncommited changes: $githaschanges
        # 
        """
        write_file(job["configfile"], job, prefix)
    end    
end


end

