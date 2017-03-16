#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re

INT_RE = "[+\-\d+]+"
FLOAT_RE = "[+-]?\d+\.\d*[eE]?[+-]?\d*"
NUMBER_RE = "[+-]?\d+\.?\d*[eE]?[+-]?\d*"

def join_re(re_iterable):
    replaced = re_iterable
    converter = re_iterable
    # substitute integer regex
    replaced = [INT_RE if item is int else item for item in replaced]
    # substitute float regex
    replaced = [FLOAT_RE if item is float else item for item in replaced]
    # add parantheses
    replaced = ["({0})".format(item) for item in replaced]
    # allow multiple whitespace between items
    joined = "\s+".join(replaced)
    joined = re.compile(joined)

    #construct converter
    converter = [item if item in (int, float) else str for item in converter]

    return joined, converter

def match_lines(lines, line_re, line_converter):
    matched = list()
    for line in lines:
        line = line.strip()
        mobj = line_re.match(line)
        if mobj:
            converted = tuple([conv(item.strip()) for item, conv in
                zip(mobj.groups(), line_converter)])
            matched.append(converted)
    return matched

def find_number(text, line_str, conv_to, limit):
    """Expects a string 'line_str'. This function then constructs
    a regular expression from this supplied string and appends the
    number regex. It matches all lines containing 'line_str' and
    returns the following floats (one per line) as a list. A limit
    is optional to control the length of the returned list."""

    line_re, converter = join_re((line_str, conv_to))
    lines = line_re.findall(text)
    return [conv_to(line[1]) for line in lines][:limit]

def find_floats(text, line_str, limit=None):
    return find_number(text, line_str, float, limit)

def find_ints(text, line_str, limit=None):
    return find_number(text, line_str, int, limit)
