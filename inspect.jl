using Gtk, Cairo, Graphics, ComfyCommons.info
using ComfyCommons.Yaml

# 
# General layout
# 

# top-level window with vertical layout
window_vbox = GtkBox(:v) # :h makes a horizontal layout, :v a vertical layout
window = GtkWindow(window_vbox, "Inspect mt-exp Plots")

# canvas
canvas = GtkCanvas(800, 320)
push!(window_vbox, canvas)
setproperty!(window_vbox, :expand, canvas, true)

# controls in a notebook, below which the full regex is presented
ctrl_vbox = GtkBox(:v)
ctrl_notebook = GtkNotebook()
push!(ctrl_vbox, ctrl_notebook)
push!(window_vbox, ctrl_vbox)
setproperty!(ctrl_vbox, :spacing, 10)
setproperty!(ctrl_vbox, :margin, 20)

# regular expression control
ctrl_regex_box = GtkBox(:h)
ctrl_regex = GtkEntry()
ctrl_regex_matches = GtkLabel("")
push!(ctrl_regex_box, ctrl_regex)
push!(ctrl_regex_box, ctrl_regex_matches)
push!(ctrl_vbox, ctrl_regex_box)
setproperty!(ctrl_regex_box, :expand, ctrl_regex, true)
setproperty!(ctrl_regex_box, :spacing, 20)
setproperty!(ctrl_regex_matches, "use-markup", true)

# copy tex-src control
ctrl_cp_box = GtkBox(:h)
ctrl_cp = GtkEntry()
ctrl_cp_button = GtkButton("Copy .tex sources")
push!(ctrl_cp_box, ctrl_cp)
push!(ctrl_cp_box, ctrl_cp_button)
push!(ctrl_vbox, ctrl_cp_box)
setproperty!(ctrl_cp_box, :expand, ctrl_cp, true)
setproperty!(ctrl_cp_box, :spacing, 20)
setproperty!(ctrl_cp, :text, "../mt/img/exp-tex/")




getreplacement(wconf::Dict{Any, Any}) =
    if wconf["type"] == "GtkEntry"
        getproperty(wconf["widget"], :text, String)
    elseif wconf["type"] == "GtkBox"
        replacements = filter(s -> !isempty(s), map(getreplacement, wconf["children"]))
        if !isempty(replacements)
            orlogic = wconf["childrenlogic"] == "or"
            joined = join(replacements, orlogic ? "|" : "")
            length(replacements) == 1 || !orlogic ? joined : "($joined)"
        else "" end
    elseif wconf["type"] == "GtkCheckButton"
        getproperty(wconf["widget"], :active, Bool) ? wconf["activetext"] : ""
    end

regexcallback(tconf::Dict{Any, Any}) = (widget) -> begin # callback building regular expression from controls
    regex = tconf["regex"]
    for wconf in tconf["widgets"]
        if contains(regex, "\$" * wconf["id"])
            regex = replace(regex, "\$" * wconf["id"], getreplacement(wconf))
        end
    end
    @async setproperty!(ctrl_regex, :text, regex) # async prevents segfault on regex_callback 
end

SIGNAL_TYPES = Dict( "GtkEntry"       => "changed",
                     "GtkCheckButton" => "clicked" )

function pushwidget!(parent::Gtk.GtkWidget, wconf::Dict{Any, Any},
                     callback::Union{Function, Void} = nothing)
    
    # constructor call
    widget = eval(parse( wconf["type"] * "(" * get(wconf, "parameters", "") * ")" ))
    wconf["widget"] = widget
    
    # add to parent
    push!(parent, widget)
    if get(wconf, "expand", false)
        setproperty!(parent, :expand, widget, true)
        setproperty!(parent, :spacing, 20)
    end
    
    # properties
    for pconf in get(wconf, "properties", [])
        setproperty!(widget, Symbol(pconf[1]), pconf[2])
    end
    
    # children (if any)
    for cconf in get(wconf, "children", [])
        pushwidget!(widget, cconf, callback)
    end
    
    # connect callback signal
    if callback != nothing && haskey(SIGNAL_TYPES, wconf["type"])
        signal_connect(callback, widget, SIGNAL_TYPES[wconf["type"]])
    end
    
end

function pushtab!(parent::Gtk.GtkWidget, tconf::Dict{Any, Any})
    
    # layout
    hbox = GtkBox(:h)
    push!(parent, hbox)
    setproperty!(hbox, :spacing, 20)
    setproperty!(hbox, :margin, 10)
    setproperty!(parent, :tab_label, hbox, tconf["label"])
    tconf["widget"] = hbox
    
    # widgets
    tcallback = regexcallback(tconf)
    for wconf in tconf["widgets"]
        pushwidget!(hbox, wconf, tcallback)
    end
    
    # tab switch callback
    signal_connect((nbleaf, leaf, index) -> if leaf == hbox tcallback(nbleaf) end,
                   parent, "switch-page")
    
end

# 
# Configured tabs
# 
conf = load_file("conf/inspect.yml")
for tconf in conf["tabs"]
    pushtab!(ctrl_notebook, tconf)
end



# 
# Canvas and Copy
# 

# update canvas in callback
function canvascallback(widget)
    
    # match files
    regex = getproperty(ctrl_regex, :text, String)
    fnames = map(f -> joinpath(dirname(regex), f), try readdir(dirname(regex)) catch String[] end)
    totalnum = length(fnames)
    fnames = filter(f -> ismatch(Regex(regex), f), fnames)
    matchlabel = "$(length(fnames)) of $totalnum files match"
    if length(fnames) > 2
        matchlabel = "<span foreground=\"red\">$matchlabel</span>"
    end
    setproperty!(ctrl_regex_matches, :label, matchlabel)
    
    # update images
    if length(fnames) > 0
        @guarded draw(canvas) do widget2
            ctx = getgc(canvas)
            w = width(canvas)
            h = height(canvas)
            Graphics.save(ctx)
            
            rectangle(ctx, 0, 0, w, h)
            set_source_rgb(ctx, 1, 1, 1)
            fill(ctx)
            
            images = read_from_png.(fnames[1:min(length(fnames), 2)])
            scaling = min(w / (2.1*maximum(map(i -> i.width, images))), h / (1.05*maximum(map(i -> i.height, images))))
            scale(ctx, scaling, scaling)
            set_source_surface(ctx, images[1], (w/scaling - 2 * images[1].width)/3, (h/scaling - images[1].height)/2)
            paint(ctx)
            
            if length(images) > 1
                set_source_surface(ctx, images[2], images[2].width + (w/scaling - 2 * images[2].width)*2/3, (h/scaling - images[2].height)/2)
                paint(ctx)
            end
            
            restore(ctx)
        end
    end
    
end
signal_connect(canvascallback, ctrl_regex, "changed")



# callback for copy button
function cp_callback(widget)
    targetdir = getproperty(ctrl_cp, :text, String)
    regex = getproperty(ctrl_regex, :text, String)
    fnames = map(f -> joinpath(dirname(regex), f), readdir(dirname(regex)))
    totalnum = length(fnames)
    fnames = filter(f -> ismatch(Regex(regex), f), fnames)
    fnames = map(f -> replace(f, ".png", ".tex"), fnames)
    fnames = map(f -> replace(f, "png/", "tex/"), fnames)
    fnames = filter(isfile, fnames)
    
    if length(fnames) > 0
        buf = IOBuffer()
        println(buf, "Do you want to copy the following files to $targetdir?")
        for f in fnames
            println(buf, "- $f")
        end
        
        # confirm in dialog
        if ask_dialog(String(take!(buf)), "Cancel", "Yes", window)
            for f in fnames
                targetpath = joinpath(targetdir, basename(f))
                println(f, " -> ", targetpath)
                cp(f, targetpath, remove_destination=true)
            end
            info_dialog(".tex sources were successfully copied!", window)
        end
    else
        buf = IOBuffer()
        println(buf, "Could not find corresponding .tex sources:")
        for f in fnames
            println(buf, "- $f")
        end
        warn_dialog(String(take!(buf)), window)
    end
end
signal_connect(cp_callback, ctrl_cp_button, "clicked")




# update window and wait for it to be closed
showall(window)
if !isinteractive()
    c = Condition()
    signal_connect(window, :destroy) do widget
        notify(c)
    end
    wait(c)
else
    nothing # do not return anything
end

