# Future Integration: lotka-volterra-agents-c

## Current State
C implementation of generalized Lotka-Volterra dynamics for multi-species agent ecology. Models predator-prey, competition, and mutualism relationships between ternary agent species.

## Integration Opportunities

### With strategy-ecology-c
Strategy ecology IS Lotka-Volterra applied to ternary strategies. The two C ports are complementary: `lotka-volterra-agents-c` provides the dynamics engine, `strategy-ecology-c` provides the strategy-specific interactions. Together: ecological dynamics + strategy competition = complete on-device population modeling.

### With ternary-cell (Population Dynamics)
Cell populations follow ecological dynamics. Cell types (explorer, stabilizer, specialist) are "species" competing for grid resources. Lotka-Volterra dynamics model the population oscillations. The C port runs on-device, predicting cell population trends for proactive GC adjustments.

### With ternary-swarm
Swarm population dynamics are ecological. PSO particle populations grow and shrink. Lotka-Volterra models the swarm size dynamics — when to spawn new particles (prey abundance), when to remove particles (predation pressure from the problem landscape).

## Potential in Mature Systems
In room-as-codespace, room populations follow ecological dynamics. Agent types compete for room access. Lotka-Volterra on the edge device predicts which agent types will thrive in a given room, enabling proactive ensign loading before demand spikes.

## Cross-Pollination Ideas
- Predator-prey as room debugging: "predator" agents hunt down malfunctioning cells
- Mutualism for room cooperation: rooms that mutually benefit should be colocated
- Competitive exclusion for resource allocation: no two rooms can occupy the same niche

## Dependencies for Next Steps
- Integration with strategy-ecology-c for combined ecological modeling
- FFI bindings for Rust ternary-cell integration
- Parameter calibration from observed room population dynamics
