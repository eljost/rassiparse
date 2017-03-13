#!/usr/bin/env python
# -*- coding: utf-8 -*-

class WithHeader:

    def as_list(self, attrs):
        return [getattr(self, attr) for attr in attrs]

    def get_headers(self, attrs):
        return [self.headers[attr] for attr in attrs]

    def as_list_with_headers(self, attrs):
        return (self.as_list(attrs), self.get_headers(attrs))
