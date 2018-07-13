#!/bin/bash
case "$1" in
    
    # 
    # execute bash shell
    # 
    bash)
        exec /bin/bash
        ;;
    
    # 
    # run julia
    # 
    julia)
        /bin/julia "${@:2}"
        ;;
    
    *)
        echo "ERROR: Missing argument in entrypoint.sh"
        exit 1
        ;;
    
esac
