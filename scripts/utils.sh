#!/bin/bash
#
# Shared utilities for SwissIsoform pipeline scripts
# Source this file in other scripts with: source utils.sh
#

# Detect if output is going to a terminal or a file
# Disable colors if output is redirected (e.g., in SLURM logs)
if [ -t 1 ]; then
    # Output is going to a terminal
    export RED='\033[0;31m'
    export GREEN='\033[0;32m'
    export YELLOW='\033[1;33m'
    export BLUE='\033[0;34m'
    export CYAN='\033[0;36m'
    export NC='\033[0m' # No Color
else
    # Output is redirected (e.g., to a file), disable colors
    export RED=''
    export GREEN=''
    export YELLOW=''
    export BLUE=''
    export CYAN=''
    export NC=''
fi

# Timing utilities
start_timer() {
    START_TIME=$(date +%s)
}

end_timer() {
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
}

format_duration() {
    local seconds=$1
    if [ $seconds -lt 60 ]; then
        echo "${seconds}s"
    elif [ $seconds -lt 3600 ]; then
        local minutes=$((seconds / 60))
        local secs=$((seconds % 60))
        echo "${minutes}m ${secs}s"
    else
        local hours=$((seconds / 3600))
        local minutes=$(((seconds % 3600) / 60))
        local secs=$((seconds % 60))
        echo "${hours}h ${minutes}m ${secs}s"
    fi
}

# Print utilities
print_header() {
    local title="$1"
    local width=64
    echo -e "${BLUE}╔$(printf '═%.0s' $(seq 1 $((width-2))))╗${NC}"
    printf "${BLUE}║${NC} %-$((width-4))s ${BLUE}║${NC}\n" "$title"
    echo -e "${BLUE}╚$(printf '═%.0s' $(seq 1 $((width-2))))╝${NC}"
}

print_section() {
    local title="$1"
    echo ""
    echo -e "${CYAN}━━━ $title ━━━${NC}"
}

print_success() {
    echo -e "${GREEN}✓${NC} $1"
}

print_error() {
    echo -e "${RED}✗${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

print_info() {
    echo -e "${YELLOW}→${NC} $1"
}

print_step() {
    echo -e "${BLUE}→${NC} $1"
}
