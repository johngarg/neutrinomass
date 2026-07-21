#!/usr/bin/env bash

# Task IDs are seed-major for the pinned 243-operator inventory. The preflight
# recomputes these expressions from the manifest and refuses to continue if
# the inventory or tier assignment has drifted.
PRIORITY4_NORMAL_TASK_ARRAY="0-167,169-172,174,179,188-189,191-410,412-415,417,422,431-432,434-485"
PRIORITY4_HIGH_TASK_ARRAY="168,173,175-178,180-187,190,411,416,418-421,423-430,433"
