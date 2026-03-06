# Dataset README

## Design (Experiment 1 and Experiment 2)
- Within-subject, repeated measures
- Factors:
  - Set size (number of bars): 3, 4, 5
  - Bar width: 0.25, 0.40, 0.55 deg
  - Bar center-to-center spacing: 0.7, 0.8, 1.1 deg
- Repetitions: 2 matched trials per condition
- Blocks: 6
- Trials per block: 54
- Total trials per participant: 324
  - `3 widths x 3 spacings x 3 set sizes x 2 repetitions x 6 blocks = 324`
- Within each block, stimulus side (left/right) was randomized, and width/spacing conditions were randomized.

## Shared Procedure
1. Participants received instructions and completed 10 practice trials.
2. Each trial:
   - Fixation point (black disc, 0.11 deg diameter) at screen center.
   - Fixation duration randomly selected from 1.00, 1.25, 1.50, 1.75, or 2.00 s.
   - Stimulus displayed for 150 ms on the left or right side of fixation.
   - Participant reported perceived number of bars using numpad keys (0-9).
   - Probe phase (see experiment-specific details below).
   - Participant pressed Enter to continue to the next trial.
3. Main experiment duration was about 45 minutes per participant.

## Experiment-Specific Probe Phase
### Experiment 1
- Probe bars reappeared at center with count equal to participant's reported number.
- Participants adjusted:
  - Bar width (step size: 0.05 deg/key press)
  - Bar spacing (step size: 0.05 deg/key press)
- Arrow-key mapping (width vs spacing control) was counterbalanced across participants.
- Initial values:
  - Width: 0.20 to 0.60 deg in 0.05 deg increments
  - Spacing: 0.30 to 1.20 deg in 0.10 deg increments

### Experiment 2
- Same as Experiment 1 except during probe adjustment:
  - Only one probe bar was presented.
  - Participants adjusted width only.
- Spacing adjustment was not performed in this experiment.


