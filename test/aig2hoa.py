#!/usr/bin/env python3
#
# Copyright Tom van Dijk
#
"""Translate an ASCII AIGER (Mealy) controller into a HOA automaton, all of
whose runs are accepting, whose language is the set of traces over the inputs
and outputs that the controller can produce. With Spot,

    aig2hoa.py controller.aag spec.ehoa > controller.hoa
    autfilt -q controller.hoa --included-in=spec.ehoa

exits with 0 iff every trace of the controller is accepted by the
specification. This is a check that shares no code with knor or aigcheck: it
has its own AIGER reader and only reads the AP: and controllable-AP: lines of
the specification. It explores the reachable latch valuations over all input
valuations, so it is only meant for controllers with few inputs.

Exits with 2 on malformed controllers (unnamed or duplicate inputs/outputs,
uninitialised latches, inputs/outputs that do not match the specification).
"""
import shlex
import sys


def die(msg):
    sys.stderr.write("aig2hoa: " + msg + "\n")
    sys.exit(2)


def read_aiger(filename):
    text = open(filename).read().split("\n")
    header = text[0].split()
    if not header or header[0] != "aag":
        die("only ASCII AIGER (aag) is supported")
    _, num_inputs, num_latches, num_outputs, num_ands = map(int, header[1:6])
    if any(int(x) for x in header[6:]):
        die("bad/constraint/justice/fairness sections are not supported")
    pos = 1
    inputs = [int(text[pos + k]) for k in range(num_inputs)]
    pos += num_inputs
    latches = []
    for k in range(num_latches):
        fields = text[pos + k].split()
        reset = int(fields[2]) if len(fields) > 2 else 0
        if reset not in (0, 1):
            die("latch %s is not initialised to 0 or 1" % fields[0])
        latches.append((int(fields[0]), int(fields[1]), reset))
    pos += num_latches
    outputs = [int(text[pos + k]) for k in range(num_outputs)]
    pos += num_outputs
    ands = [tuple(map(int, text[pos + k].split())) for k in range(num_ands)]
    pos += num_ands
    input_names, output_names = [None] * num_inputs, [None] * num_outputs
    for line in text[pos:]:
        if line.startswith("c"):
            break
        if line.startswith("i"):
            k, name = line[1:].split(" ", 1)
            input_names[int(k)] = name
        elif line.startswith("o"):
            k, name = line[1:].split(" ", 1)
            output_names[int(k)] = name
    if None in input_names or None in output_names:
        die("every input and output must be named")
    if (len(set(input_names)) != num_inputs or len(set(output_names)) != num_outputs
            or set(input_names) & set(output_names)):
        die("the input and output names must be distinct")
    return inputs, latches, outputs, ands, input_names, output_names


def read_spec_aps(filename):
    aps, controllable = None, None
    for line in open(filename):
        if line.startswith("AP:"):
            aps = shlex.split(line)[2:]
        elif line.startswith("controllable-AP:"):
            controllable = [int(x) for x in line.split()[1:]]
        elif line.startswith("--BODY--"):
            break
    if aps is None or controllable is None:
        die("the specification has no AP: or controllable-AP: line")
    return aps, controllable


def main():
    if len(sys.argv) != 3:
        sys.stderr.write("usage: aig2hoa.py controller.aag spec.ehoa\n")
        sys.exit(2)
    inputs, latches, outputs, ands, input_names, output_names = read_aiger(sys.argv[1])
    aps, controllable = read_spec_aps(sys.argv[2])
    if sorted(input_names) != sorted(aps[k] for k in range(len(aps)) if k not in controllable):
        die("the inputs %s are not the uncontrollable APs" % input_names)
    if sorted(output_names) != sorted(aps[k] for k in controllable):
        die("the outputs %s are not the controllable APs" % output_names)
    ap_index = {name: k for k, name in enumerate(aps)}

    # evaluate the AND gates in an order where their arguments are known (ASCII AIGER allows any order)
    definition = {lhs >> 1: (lhs, rhs0, rhs1) for lhs, rhs0, rhs1 in ands}
    done = {0} | {lit >> 1 for lit in inputs} | {lit >> 1 for lit, _, _ in latches}
    order, active = [], set()
    for root in definition:
        stack = [(root, False)]
        while stack:
            var, expanded = stack.pop()
            if expanded:
                active.discard(var)
                done.add(var)
                order.append(definition[var])
                continue
            if var in done:
                continue
            if var in active:
                die("cyclic AND gates")
            if var not in definition:
                die("undefined literal %d" % (2 * var))
            active.add(var)
            stack.append((var, True))
            for arg in definition[var][1:]:
                if (arg >> 1) not in done:
                    stack.append((arg >> 1, False))

    def step(input_values, latch_values):
        value = {0: 0}
        for lit, b in zip(inputs, input_values):
            value[lit >> 1] = b
        for (lit, _, _), b in zip(latches, latch_values):
            value[lit >> 1] = b
        lit_value = lambda lit: value[lit >> 1] ^ (lit & 1)
        for lhs, rhs0, rhs1 in order:
            value[lhs >> 1] = lit_value(rhs0) & lit_value(rhs1)
        return [lit_value(o) for o in outputs], tuple(lit_value(nxt) for _, nxt, _ in latches)

    initial = tuple(reset for _, _, reset in latches)
    states = {initial: 0}
    queue = [initial]
    edges = []   # per state: (label, successor)
    for latch_values in queue:
        out = []
        for u in range(1 << len(inputs)):
            input_values = [(u >> b) & 1 for b in range(len(inputs))]
            output_values, successor = step(input_values, latch_values)
            if successor not in states:
                states[successor] = len(queue)
                queue.append(successor)
            literals = ["%s%d" % ("" if b else "!", ap_index[name]) for name, b in zip(input_names, input_values)]
            literals += ["%s%d" % ("" if b else "!", ap_index[name]) for name, b in zip(output_names, output_values)]
            out.append(("&".join(literals) if literals else "t", states[successor]))
        edges.append(out)

    w = sys.stdout.write
    w("HOA: v1\nStates: %d\nStart: 0\n" % len(queue))
    w("AP: %d %s\n" % (len(aps), " ".join('"%s"' % a for a in aps)))
    w("acc-name: all\nAcceptance: 0 t\nproperties: trans-labels explicit-labels deterministic\n--BODY--\n")
    for s, out in enumerate(edges):
        w("State: %d\n" % s)
        for label, successor in out:
            w("[%s] %d\n" % (label, successor))
    w("--END--\n")


if __name__ == "__main__":
    main()
