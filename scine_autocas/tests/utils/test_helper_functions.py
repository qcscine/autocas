# -*- coding: utf-8 -*-
__copyright__ = """ This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import unittest
from dataclasses import dataclass

from scine_autocas.utils import helper_functions


@dataclass
class TestStaticClass:
    __test__ = False
    x = "x"

    @dataclass
    class TestStaticSubClass:
        __test__ = False
        sub = "sub"

    @dataclass
    class TestOtherStaticSubClass:
        __test__ = False
        sub2 = "sub2"


class TestSubSubClass:
    __test__ = False

    def __init__(self):
        self.sub_sub_class_attr = "subsub"


class TestSubClass:
    __test__ = False

    def __init__(self):
        self.sub_class_attr = "sub"
        self.sub_sub_class = TestSubSubClass()


class TestClass:
    __test__ = False

    def __init__(self, test_attr: str = "haha"):
        self.a = 5
        self.b = "r"
        self.c = [5, 4, 6]
        self.d = {"abc": 5}
        self.e = ("xyz", "def")
        self.test_attr = test_attr
        self.sub_class = TestSubClass()


class TestHelperFunctions(unittest.TestCase):
    def setUp(self):
        self.yaml_name = "test.yaml"

        self.test_dict = {
            "TestClass": {
                "a": 5,
                "b": "r",
                "c": [5, 4, 6],
                "d": {"abc": 5},
                "e": ("xyz", "def"),
                "test_attr": "haha",
                "sub_class": {
                    "TestSubClass": {
                        "sub_class_attr": "sub",
                        "sub_sub_class": {
                            "TestSubSubClass": {
                                "sub_sub_class_attr": "subsub"
                            }
                        }
                    }
                }
            }
        }

        self.test_static_dict = {
            "TestStaticClass": {
                "x": "x",
                "TestStaticSubClass": {
                    "sub": "sub"
                },
                "TestOtherStaticSubClass": {
                    "sub2": "sub2"
                }
            }
        }

    def tearDown(self):
        try:
            os.remove(self.yaml_name)
        except FileNotFoundError:
            pass

    def test_get_parameters(self):
        test_dict = helper_functions.get_all_params(TestStaticClass)
        self.assertEqual(self.test_static_dict, test_dict)

        TestStaticClass.TestStaticSubClass.sub = "bus"
        test_dict = helper_functions.get_all_params(TestStaticClass)
        self.test_static_dict["TestStaticClass"]["TestStaticSubClass"]["sub"] = "bus"
        self.assertEqual(self.test_static_dict, test_dict)

        test_obj = TestClass()
        test_dict = helper_functions.get_all_params(test_obj)
        self.assertEqual(self.test_dict, test_dict)

        test_obj = TestClass(test_attr="jkl")
        test_dict = helper_functions.get_all_params(test_obj)
        self.test_dict["TestClass"]["test_attr"] = "jkl"
        self.assertEqual(self.test_dict, test_dict)

    def test_set_parameters(self):
        self.test_static_dict["TestStaticClass"]["TestStaticSubClass"]["sub"] = "bus"
        helper_functions.set_all_params(TestStaticClass, self.test_static_dict)
        test_dict = helper_functions.get_all_params(TestStaticClass)
        self.assertEqual(self.test_static_dict, test_dict)

        test_obj = TestClass()
        self.test_dict["TestClass"]["test_attr"] = "jkl"
        test_obj = helper_functions.set_all_params(test_obj, self.test_dict)
        test_dict = helper_functions.get_all_params(test_obj)
        self.assertEqual(self.test_dict, test_dict)


if __name__ == "__main__":
    unittest.main()
