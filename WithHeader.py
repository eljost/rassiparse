#!/usr/bin/env python
# -*- coding: utf-8 -*-

class WithHeader:

    hartree2eV = 27.211386
    hartree2nm = 45.5640

    def as_list(self, attrs):
        return [getattr(self, attr) for attr in attrs]

    def as_str_list(self, attrs):
        fmts = self.get_formats(attrs)
        return [fmt.format(getattr(self, attr)) for attr, fmt
                in zip(attrs, fmts)]

    def get_headers(self, attrs):
        return [self.headers[attr][0] for attr in attrs]

    def get_formats(self, attrs):
        return [self.headers[attr][1] for attr in attrs]

    def as_list_with_headers(self, attrs):
        return (self.as_list(attrs), self.get_headers(attrs))
