#!/bin/bash

# Source the bash configuration to load any environment variables or aliases
. ~/.bashrc

# Default behavior: Run an interactive bash session if no arguments are provided
if [ -z "$*" ]; then
    exec /bin/bash
elif [ "$1" == "trimmomatic" ]; then
    # Shift the command-line arguments and execute Trimmomatic
    shift
    exec java -jar $DEST/$VERSION/$SW_NAME-$VERSION/trimmomatic-0.39.jar "$@"
elif [ "$1" == "STAR" ]; then
    # Shift the command-line arguments and execute STAR
    shift
    exec /home/bin/STAR "$@"
else
    # For any other command, execute it as provided
    exec "$@"
fi
