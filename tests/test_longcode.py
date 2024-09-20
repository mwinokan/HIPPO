#!/usr/bin/env python

import os
import sys

sys.path.insert(0, os.path.abspath("../"))

import unittest

from hippo.fragalysis import (
    parse_observation_longcode,
    find_observation_longcode_matches,
)

from mlog import setup_logger

logger = setup_logger("TestLongCode")

codes = [
    "CHIKV_MacB-x0270_A_304_1_CHIKV_MacB-x2352+C+401+1",
    "CHIKV_MacB-x0289_C_401_CHIKV_MacB-x0300+A+401+1",
    "CHIKV_MacB-x0300_A_401_CHIKV_MacB-x0300+A+401+1",
    "6vuq_A_500_1_CHIKV_MacB-x2352+C+401+1",
    "6w8q_A_201_1_6w8q+A+201+1",
    "CHIKV_MacB-x0270_A_304_1_CHIKV_MacB-x2352+C+401+1",
    "CHIKV_MacB-x0406_A_401_CHIKV_MacB-x0300+A+401+1",
]


class TestLongCode(unittest.TestCase):

    def test_longcodes(self):
        logger.header("test_longcodes")

        for code in codes:
            result = parse_observation_longcode(code)

            print(result)

    def test_match(self):
        logger.header("test_match")

        for code in codes:
            matches = find_observation_longcode_matches(code, list(EXAMPLE_MAP.keys()))

            logger.var(code, matches)

            if len(matches) != 1:
                raise Exception(code)


EXAMPLE_MAP = {
    "6vuq_A_500_1_CHIKV_MacB-x2352+C+401+1": 6221,
    "6vuq_B_200_1_CHIKV_MacB-x2352+C+401+1": 6222,
    "6vuq_C_200_1_CHIKV_MacB-x2352+C+401+1": 6223,
    "6vuq_D_200_1_CHIKV_MacB-x2352+C+401+1": 6224,
    "6w0t_A_201_1_CHIKV_MacB-x2352+C+401+1": 6201,
    "6w0t_B_201_1_CHIKV_MacB-x2352+C+401+1": 6202,
    "6w0t_C_201_1_CHIKV_MacB-x2352+C+401+1": 6203,
    "6w0t_D_201_1_CHIKV_MacB-x2352+C+401+1": 6204,
    "6w7h_A_201_1_CHIKV_MacB-x2352+C+401+1": 6217,
    "6w7h_B_201_1_CHIKV_MacB-x2352+C+401+1": 6218,
    "6w7h_C_201_1_CHIKV_MacB-x2352+C+401+1": 6219,
    "6w7h_D_201_1_CHIKV_MacB-x2352+C+401+1": 6220,
    "6w8k_A_201_1_CHIKV_MacB-x2352+C+401+1": 6205,
    "6w8k_B_201_1_CHIKV_MacB-x2352+C+401+1": 6206,
    "6w8k_C_201_1_CHIKV_MacB-x2352+C+401+1": 6207,
    "6w8k_D_201_1_CHIKV_MacB-x2352+C+401+1": 6208,
    "6w8m_A_201_1_CHIKV_MacB-x2352+C+401+1": 6209,
    "6w8m_B_201_1_CHIKV_MacB-x2352+C+401+1": 6210,
    "6w8m_C_201_1_CHIKV_MacB-x2352+C+401+1": 6211,
    "6w8m_D_201_1_CHIKV_MacB-x2352+C+401+1": 6212,
    "6w8q_A_201_1_6w8q+A+201+1": 5822,
    "6w8q_B_201_1_CHIKV_MacB-x2352+C+401+1": 5819,
    "6w8q_C_201_1_CHIKV_MacB-x2352+C+401+1": 5820,
    "6w8q_D_201_1_CHIKV_MacB-x2352+C+401+1": 5821,
    "6w8y_A_201_1_CHIKV_MacB-x2352+C+401+1": 6214,
    "6w8y_B_201_1_CHIKV_MacB-x2352+C+401+1": 6215,
    "6w8y_C_201_1_CHIKV_MacB-x2352+C+401+1": 6216,
    "6w8y_D_201_1_CHIKV_MacB-x2352+C+401+1": 6213,
    "6w8z_A_201_1_CHIKV_MacB-x2352+C+401+1": 6193,
    "6w8z_B_201_1_CHIKV_MacB-x2352+C+401+1": 6194,
    "6w8z_C_201_1_CHIKV_MacB-x2352+C+401+1": 6195,
    "6w8z_D_201_1_CHIKV_MacB-x2352+C+401+1": 6196,
    "6w91_A_201_1_CHIKV_MacB-x2352+C+401+1": 6197,
    "6w91_B_201_1_CHIKV_MacB-x2352+C+401+1": 6198,
    "6w91_C_201_1_CHIKV_MacB-x2352+C+401+1": 6199,
    "6w91_D_201_1_CHIKV_MacB-x2352+C+401+1": 6200,
    "CHIKV_MacB-x0270_A_304_1_CHIKV_MacB-x2352+C+401+1": 5802,
    "CHIKV_MacB-x0281_C_304_1_CHIKV_MacB-x2352+C+401+1": 5803,
    "CHIKV_MacB-x0289_C_401_1_CHIKV_MacB-x2352+C+401+1": 5804,
    "CHIKV_MacB-x0289_C_401_CHIKV_MacB-x0300+A+401+1": 34,
    "CHIKV_MacB-x0294_A_401_1_CHIKV_MacB-x2352+C+401+1": 5806,
    "CHIKV_MacB-x0294_A_501_1_CHIKV_MacB-x0294+A+501+1": 5805,
    "CHIKV_MacB-x0294_A_501_CHIKV_MacB-x0294+A+501+1": 63,
    "CHIKV_MacB-x0294_B_304_1_CHIKV_MacB-x2352+C+401+1": 5807,
    "CHIKV_MacB-x0294_B_304_CHIKV_MacB-x0300+A+401+1": 47,
    "CHIKV_MacB-x0295_A_401_1_CHIKV_MacB-x2352+C+401+1": 5808,
    "CHIKV_MacB-x0295_A_401_CHIKV_MacB-x0300+A+401+1": 72,
    "CHIKV_MacB-x0300_A_401_1_CHIKV_MacB-x2352+C+401+1": 5809,
    "CHIKV_MacB-x0300_C_304_1_CHIKV_MacB-x2352+C+401+1": 5810,
    "CHIKV_MacB-x0300_C_304_CHIKV_MacB-x0300+A+401+1": 52,
    "CHIKV_MacB-x0305_C_304_1_CHIKV_MacB-x2352+C+401+1": 5811,
    "CHIKV_MacB-x0305_C_304_CHIKV_MacB-x0300+A+401+1": 68,
    "CHIKV_MacB-x0312_A_304_1_CHIKV_MacB-x2352+C+401+1": 5814,
    "CHIKV_MacB-x0312_C_401_1_CHIKV_MacB-x2352+C+401+1": 5812,
    "CHIKV_MacB-x0312_C_401_CHIKV_MacB-x0300+A+401+1": 24,
    "CHIKV_MacB-x0312_D_401_1_CHIKV_MacB-x2352+C+401+1": 5813,
    "CHIKV_MacB-x0312_D_401_CHIKV_MacB-x0300+A+401+1": 2,
    "CHIKV_MacB-x0314_B_401_1_CHIKV_MacB-x2352+C+401+1": 5817,
    "CHIKV_MacB-x0314_B_401_CHIKV_MacB-x0300+A+401+1": 9,
    "CHIKV_MacB-x0314_C_304_1_CHIKV_MacB-x2352+C+401+1": 5815,
    "CHIKV_MacB-x0314_D_401_1_CHIKV_MacB-x2352+C+401+1": 5816,
    "CHIKV_MacB-x0316_C_304_1_CHIKV_MacB-x2352+C+401+1": 5818,
    "CHIKV_MacB-x0317_C_401_1_CHIKV_MacB-x2352+C+401+1": 5823,
    "CHIKV_MacB-x0317_C_401_CHIKV_MacB-x0300+A+401+1": 28,
    "CHIKV_MacB-x0353_A_304_1_CHIKV_MacB-x2352+C+401+1": 5824,
    "CHIKV_MacB-x0353_B_304_1_CHIKV_MacB-x2352+C+401+1": 5825,
    "CHIKV_MacB-x0394_D_304_1_CHIKV_MacB-x2352+C+401+1": 5826,
    "CHIKV_MacB-x0398_C_401_1_CHIKV_MacB-x2352+C+401+1": 5827,
    "CHIKV_MacB-x0398_C_401_CHIKV_MacB-x0300+A+401+1": 27,
    "CHIKV_MacB-x0399_A_401_1_CHIKV_MacB-x2352+C+401+1": 5828,
    "CHIKV_MacB-x0399_A_402_1_CHIKV_MacB-x2352+C+401+1": 5829,
    "CHIKV_MacB-x0399_B_401_1_CHIKV_MacB-x2352+C+401+1": 5830,
    "CHIKV_MacB-x0399_B_402_1_CHIKV_MacB-x2352+C+401+1": 5831,
    "CHIKV_MacB-x0399_C_308_1_CHIKV_MacB-x2352+C+401+1": 5833,
    "CHIKV_MacB-x0399_C_309_1_CHIKV_MacB-x2352+C+401+1": 5832,
    "CHIKV_MacB-x0399_D_304_1_CHIKV_MacB-x2352+C+401+1": 5834,
    "CHIKV_MacB-x0404_A_401_1_CHIKV_MacB-x2352+C+401+1": 5835,
    "CHIKV_MacB-x0404_A_401_CHIKV_MacB-x0300+A+401+1": 65,
    "CHIKV_MacB-x0404_C_304_1_CHIKV_MacB-x2352+C+401+1": 5836,
    "CHIKV_MacB-x0406_A_401_1_CHIKV_MacB-x2352+C+401+1": 5837,
    "CHIKV_MacB-x0406_B_304_1_CHIKV_MacB-x2352+C+401+1": 5838,
    "CHIKV_MacB-x0436_A_401_1_CHIKV_MacB-x2352+C+401+1": 5839,
    "CHIKV_MacB-x0436_C_304_1_CHIKV_MacB-x2352+C+401+1": 5840,
    "CHIKV_MacB-x0436_C_304_CHIKV_MacB-x0300+A+401+1": 29,
    "CHIKV_MacB-x0441_D_401_1_CHIKV_MacB-x2352+C+401+1": 5841,
    "CHIKV_MacB-x0477_A_401_1_CHIKV_MacB-x2352+C+401+1": 5842,
    "CHIKV_MacB-x0477_D_401_1_CHIKV_MacB-x2352+C+401+1": 5843,
    "CHIKV_MacB-x0505_A_401_1_CHIKV_MacB-x2352+C+401+1": 5845,
    "CHIKV_MacB-x0505_B_401_1_CHIKV_MacB-x2352+C+401+1": 5846,
    "CHIKV_MacB-x0505_C_401_1_CHIKV_MacB-x2352+C+401+1": 5847,
    "CHIKV_MacB-x0522_A_401_1_CHIKV_MacB-x2352+C+401+1": 5848,
    "CHIKV_MacB-x0522_A_401_CHIKV_MacB-x0300+A+401+1": 23,
    "CHIKV_MacB-x0522_B_401_1_CHIKV_MacB-x2352+C+401+1": 5849,
    "CHIKV_MacB-x0522_D_306_1_CHIKV_MacB-x2352+C+401+1": 5850,
    "CHIKV_MacB-x0540_C_401_1_CHIKV_MacB-x2352+C+401+1": 5851,
    "CHIKV_MacB-x0543_A_304_1_CHIKV_MacB-x2352+C+401+1": 5853,
    "CHIKV_MacB-x0543_B_308_1_CHIKV_MacB-x0543+B+308+1": 5852,
    "CHIKV_MacB-x0543_C_401_1_CHIKV_MacB-x2352+C+401+1": 5854,
    "CHIKV_MacB-x0544_D_401_1_CHIKV_MacB-x2352+C+401+1": 5855,
    "CHIKV_MacB-x0544_D_401_CHIKV_MacB-x0300+A+401+1": 16,
    "CHIKV_MacB-x0553_A_401_1_CHIKV_MacB-x2352+C+401+1": 5856,
    "CHIKV_MacB-x0553_B_304_1_CHIKV_MacB-x2352+C+401+1": 5857,
    "CHIKV_MacB-x0553_B_304_CHIKV_MacB-x0300+A+401+1": 59,
    "CHIKV_MacB-x0585_B_308_1_CHIKV_MacB-x2352+C+401+1": 5858,
    "CHIKV_MacB-x0585_C_401_1_CHIKV_MacB-x2352+C+401+1": 5859,
    "CHIKV_MacB-x0585_D_306_1_CHIKV_MacB-x2352+C+401+1": 5860,
    "CHIKV_MacB-x0620_C_401_1_CHIKV_MacB-x2352+C+401+1": 5861,
    "CHIKV_MacB-x0620_C_401_CHIKV_MacB-x0300+A+401+1": 22,
    "CHIKV_MacB-x0632_A_401_1_CHIKV_MacB-x2352+C+401+1": 5862,
    "CHIKV_MacB-x0632_A_401_CHIKV_MacB-x0300+A+401+1": 49,
    "CHIKV_MacB-x0692_C_401_1_CHIKV_MacB-x2352+C+401+1": 5864,
    "CHIKV_MacB-x0692_C_401_CHIKV_MacB-x0300+A+401+1": 54,
    "CHIKV_MacB-x0692_D_304_1_CHIKV_MacB-x0692+D+304+1": 5863,
    "CHIKV_MacB-x0721_B_401_1_CHIKV_MacB-x2352+C+401+1": 5865,
    "CHIKV_MacB-x0732_C_401_1_CHIKV_MacB-x2352+C+401+1": 5866,
    "CHIKV_MacB-x0732_C_401_CHIKV_MacB-x0300+A+401+1": 5,
    "CHIKV_MacB-x0756_A_304_1_CHIKV_MacB-x2352+C+401+1": 5868,
    "CHIKV_MacB-x0756_C_401_1_CHIKV_MacB-x0543+B+308+1": 5867,
    "CHIKV_MacB-x0756_C_402_1_CHIKV_MacB-x2352+C+401+1": 5869,
    "CHIKV_MacB-x0756_D_304_1_CHIKV_MacB-x2352+C+401+1": 5870,
    "CHIKV_MacB-x0760_A_401_1_CHIKV_MacB-x2352+C+401+1": 5871,
    "CHIKV_MacB-x0760_A_401_CHIKV_MacB-x0300+A+401+1": 60,
    "CHIKV_MacB-x0760_B_304_1_CHIKV_MacB-x2352+C+401+1": 5872,
    "CHIKV_MacB-x0760_C_308_1_CHIKV_MacB-x2352+C+401+1": 5873,
    "CHIKV_MacB-x0760_D_401_1_CHIKV_MacB-x2352+C+401+1": 5874,
    "CHIKV_MacB-x0822_A_401_1_CHIKV_MacB-x2352+C+401+1": 5875,
    "CHIKV_MacB-x0824_A_401_1_CHIKV_MacB-x2352+C+401+1": 5876,
    "CHIKV_MacB-x0824_A_401_CHIKV_MacB-x0300+A+401+1": 55,
    "CHIKV_MacB-x0839_A_401_1_CHIKV_MacB-x2352+C+401+1": 5877,
    "CHIKV_MacB-x0839_B_401_1_CHIKV_MacB-x2352+C+401+1": 5878,
    "CHIKV_MacB-x0839_B_401_CHIKV_MacB-x0300+A+401+1": 17,
    "CHIKV_MacB-x0839_C_308_1_CHIKV_MacB-x2352+C+401+1": 5879,
    "CHIKV_MacB-x0839_C_308_CHIKV_MacB-x0300+A+401+1": 7,
    "CHIKV_MacB-x0839_D_304_1_CHIKV_MacB-x2352+C+401+1": 5881,
    "CHIKV_MacB-x0839_D_304_CHIKV_MacB-x0300+A+401+1": 8,
    "CHIKV_MacB-x0839_D_401_1_CHIKV_MacB-x2352+C+401+1": 5880,
    "CHIKV_MacB-x0839_D_401_CHIKV_MacB-x0300+A+401+1": 15,
    "CHIKV_MacB-x0858_A_401_1_CHIKV_MacB-x2352+C+401+1": 5884,
    "CHIKV_MacB-x0858_B_304_1_CHIKV_MacB-x2352+C+401+1": 5885,
    "CHIKV_MacB-x0858_C_308_1_CHIKV_MacB-x2352+C+401+1": 5886,
    "CHIKV_MacB-x0858_D_306_1_CHIKV_MacB-x2352+C+401+1": 5887,
    "CHIKV_MacB-x0864_A_401_1_CHIKV_MacB-x2352+C+401+1": 5891,
    "CHIKV_MacB-x0864_A_401_CHIKV_MacB-x0300+A+401+1": 43,
    "CHIKV_MacB-x0864_A_501_1_CHIKV_MacB-x2352+C+401+1": 5892,
    "CHIKV_MacB-x0864_C_304_1_CHIKV_MacB-x2352+C+401+1": 5888,
    "CHIKV_MacB-x0864_C_401_1_CHIKV_MacB-x2352+C+401+1": 5889,
    "CHIKV_MacB-x0864_D_401_1_CHIKV_MacB-x2352+C+401+1": 5890,
    "CHIKV_MacB-x0888_A_304_1_CHIKV_MacB-x0543+B+308+1": 5893,
    "CHIKV_MacB-x0888_C_401_1_CHIKV_MacB-x2352+C+401+1": 5894,
    "CHIKV_MacB-x0892_A_304_1_CHIKV_MacB-x2352+C+401+1": 5895,
    "CHIKV_MacB-x0892_D_306_1_CHIKV_MacB-x2352+C+401+1": 5896,
    "CHIKV_MacB-x0894_A_401_1_CHIKV_MacB-x2352+C+401+1": 5897,
    "CHIKV_MacB-x0894_A_402_1_CHIKV_MacB-x2352+C+401+1": 5898,
    "CHIKV_MacB-x0894_B_308_1_CHIKV_MacB-x2352+C+401+1": 5901,
    "CHIKV_MacB-x0894_B_309_1_CHIKV_MacB-x2352+C+401+1": 5899,
    "CHIKV_MacB-x0894_C_308_1_CHIKV_MacB-x2352+C+401+1": 5900,
    "CHIKV_MacB-x0894_C_309_1_CHIKV_MacB-x2352+C+401+1": 5902,
    "CHIKV_MacB-x0926_A_401_1_CHIKV_MacB-x2352+C+401+1": 5903,
    "CHIKV_MacB-x0926_B_308_1_CHIKV_MacB-x2352+C+401+1": 5904,
    "CHIKV_MacB-x0926_D_304_1_CHIKV_MacB-x2352+C+401+1": 5905,
    "CHIKV_MacB-x0927_B_307_1_CHIKV_MacB-x2352+C+401+1": 5906,
    "CHIKV_MacB-x0927_C_401_1_CHIKV_MacB-x2352+C+401+1": 5907,
    "CHIKV_MacB-x0927_D_401_1_CHIKV_MacB-x2352+C+401+1": 5908,
    "CHIKV_MacB-x0934_A_304_1_CHIKV_MacB-x2352+C+401+1": 5909,
    "CHIKV_MacB-x0934_A_401_1_CHIKV_MacB-x2352+C+401+1": 5910,
    "CHIKV_MacB-x0934_C_308_1_CHIKV_MacB-x2352+C+401+1": 5911,
    "CHIKV_MacB-x0934_D_306_1_CHIKV_MacB-x2352+C+401+1": 5912,
    "CHIKV_MacB-x0935_A_401_1_CHIKV_MacB-x0692+D+304+1": 5913,
    "CHIKV_MacB-x0935_B_308_1_CHIKV_MacB-x0692+D+304+1": 5914,
    "CHIKV_MacB-x0935_C_401_1_CHIKV_MacB-x2352+C+401+1": 5916,
    "CHIKV_MacB-x0935_C_401_CHIKV_MacB-x0300+A+401+1": 50,
    "CHIKV_MacB-x0935_D_304_1_CHIKV_MacB-x0692+D+304+1": 5915,
    "CHIKV_MacB-x0969_A_304_1_CHIKV_MacB-x2352+C+401+1": 5795,
    "CHIKV_MacB-x0969_A_401_1_CHIKV_MacB-x2352+C+401+1": 5791,
    "CHIKV_MacB-x0969_C_401_1_CHIKV_MacB-x2352+C+401+1": 5792,
    "CHIKV_MacB-x0969_C_501_1_CHIKV_MacB-x2352+C+401+1": 5793,
    "CHIKV_MacB-x0969_D_401_1_CHIKV_MacB-x2352+C+401+1": 5796,
    "CHIKV_MacB-x0969_D_501_1_CHIKV_MacB-x2352+C+401+1": 5794,
    "CHIKV_MacB-x0997_B_308_1_CHIKV_MacB-x2352+C+401+1": 5789,
    "CHIKV_MacB-x0997_C_401_1_CHIKV_MacB-x2352+C+401+1": 5790,
    "CHIKV_MacB-x1075_A_304_1_CHIKV_MacB-x2352+C+401+1": 5799,
    "CHIKV_MacB-x1075_B_304_1_CHIKV_MacB-x2352+C+401+1": 5800,
    "CHIKV_MacB-x1075_D_401_1_CHIKV_MacB-x2352+C+401+1": 5801,
    "CHIKV_MacB-x1075_D_401_CHIKV_MacB-x0300+A+401+1": 26,
    "CHIKV_MacB-x1076_B_401_1_CHIKV_MacB-x2352+C+401+1": 5844,
    "CHIKV_MacB-x1091_C_308_1_CHIKV_MacB-x2352+C+401+1": 5882,
    "CHIKV_MacB-x1091_D_401_1_CHIKV_MacB-x2352+C+401+1": 5883,
    "CHIKV_MacB-x1103_B_401_1_CHIKV_MacB-x2352+C+401+1": 5919,
    "CHIKV_MacB-x1103_C_308_1_CHIKV_MacB-x2352+C+401+1": 5917,
    "CHIKV_MacB-x1103_D_304_1_CHIKV_MacB-x2352+C+401+1": 5918,
    "CHIKV_MacB-x1108_A_304_1_CHIKV_MacB-x2352+C+401+1": 5920,
    "CHIKV_MacB-x1108_B_308_1_CHIKV_MacB-x2352+C+401+1": 5921,
    "CHIKV_MacB-x1108_B_308_CHIKV_MacB-x0300+A+401+1": 3,
    "CHIKV_MacB-x1108_B_309_1_CHIKV_MacB-x2352+C+401+1": 5922,
    "CHIKV_MacB-x1108_B_309_CHIKV_MacB-x0300+A+401+1": 25,
    "CHIKV_MacB-x1108_C_308_1_CHIKV_MacB-x2352+C+401+1": 5924,
    "CHIKV_MacB-x1108_C_308_CHIKV_MacB-x0300+A+401+1": 20,
    "CHIKV_MacB-x1108_D_306_1_CHIKV_MacB-x2352+C+401+1": 5923,
    "CHIKV_MacB-x1108_D_306_CHIKV_MacB-x0300+A+401+1": 6,
    "CHIKV_MacB-x1112_A_305_1_CHIKV_MacB-x2352+C+401+1": 5925,
    "CHIKV_MacB-x1112_A_305_CHIKV_MacB-x0300+A+401+1": 74,
    "CHIKV_MacB-x1112_A_306_1_CHIKV_MacB-x2352+C+401+1": 5926,
    "CHIKV_MacB-x1112_B_304_1_CHIKV_MacB-x2352+C+401+1": 5930,
    "CHIKV_MacB-x1112_B_401_1_CHIKV_MacB-x2352+C+401+1": 5927,
    "CHIKV_MacB-x1112_C_401_1_CHIKV_MacB-x2352+C+401+1": 5931,
    "CHIKV_MacB-x1112_C_402_1_CHIKV_MacB-x2352+C+401+1": 5928,
    "CHIKV_MacB-x1112_C_402_CHIKV_MacB-x0300+A+401+1": 71,
    "CHIKV_MacB-x1112_D_306_1_CHIKV_MacB-x2352+C+401+1": 5929,
    "CHIKV_MacB-x1114_A_401_1_CHIKV_MacB-x2352+C+401+1": 5936,
    "CHIKV_MacB-x1114_B_304_1_CHIKV_MacB-x2352+C+401+1": 5932,
    "CHIKV_MacB-x1114_B_308_1_CHIKV_MacB-x2352+C+401+1": 5933,
    "CHIKV_MacB-x1114_B_308_CHIKV_MacB-x0300+A+401+1": 79,
    "CHIKV_MacB-x1114_C_308_1_CHIKV_MacB-x2352+C+401+1": 5934,
    "CHIKV_MacB-x1114_C_308_CHIKV_MacB-x0300+A+401+1": 78,
    "CHIKV_MacB-x1114_D_401_1_CHIKV_MacB-x2352+C+401+1": 5935,
    "CHIKV_MacB-x1116_A_304_1_CHIKV_MacB-x2352+C+401+1": 5940,
    "CHIKV_MacB-x1116_A_304_CHIKV_MacB-x0300+A+401+1": 56,
    "CHIKV_MacB-x1116_B_401_1_CHIKV_MacB-x2352+C+401+1": 5941,
    "CHIKV_MacB-x1116_B_402_1_CHIKV_MacB-x2352+C+401+1": 5937,
    "CHIKV_MacB-x1116_D_304_1_CHIKV_MacB-x2352+C+401+1": 5938,
    "CHIKV_MacB-x1116_D_401_1_CHIKV_MacB-x2352+C+401+1": 5939,
    "CHIKV_MacB-x1118_A_304_1_CHIKV_MacB-x2352+C+401+1": 5949,
    "CHIKV_MacB-x1118_A_304_CHIKV_MacB-x0300+A+401+1": 83,
    "CHIKV_MacB-x1118_A_305_1_CHIKV_MacB-x2352+C+401+1": 5942,
    "CHIKV_MacB-x1118_B_401_1_CHIKV_MacB-x2352+C+401+1": 5943,
    "CHIKV_MacB-x1118_B_402_1_CHIKV_MacB-x2352+C+401+1": 5950,
    "CHIKV_MacB-x1118_C_307_1_CHIKV_MacB-x2352+C+401+1": 5944,
    "CHIKV_MacB-x1118_C_308_1_CHIKV_MacB-x0692+D+304+1": 5947,
    "CHIKV_MacB-x1118_C_308_CHIKV_MacB-x0692+D+304+1": 87,
    "CHIKV_MacB-x1118_D_308_1_CHIKV_MacB-x2352+C+401+1": 5951,
    "CHIKV_MacB-x1118_D_401_1_CHIKV_MacB-x2352+C+401+1": 5945,
    "CHIKV_MacB-x1118_D_401_CHIKV_MacB-x0300+A+401+1": 80,
    "CHIKV_MacB-x1118_D_402_1_CHIKV_MacB-x2352+C+401+1": 5946,
    "CHIKV_MacB-x1118_D_403_1_CHIKV_MacB-x0692+D+304+1": 5948,
    "CHIKV_MacB-x1123_A_304_1_CHIKV_MacB-x2352+C+401+1": 5958,
    "CHIKV_MacB-x1123_A_305_1_CHIKV_MacB-x2352+C+401+1": 5959,
    "CHIKV_MacB-x1123_B_304_1_CHIKV_MacB-x2352+C+401+1": 5960,
    "CHIKV_MacB-x1123_B_308_1_CHIKV_MacB-x2352+C+401+1": 5967,
    "CHIKV_MacB-x1123_B_309_1_CHIKV_MacB-x2352+C+401+1": 5968,
    "CHIKV_MacB-x1123_B_310_1_CHIKV_MacB-x1123+B+310+1": 5966,
    "CHIKV_MacB-x1123_C_401_1_6w8q+A+201+1": 5964,
    "CHIKV_MacB-x1123_C_402_1_CHIKV_MacB-x2352+C+401+1": 5961,
    "CHIKV_MacB-x1123_C_403_1_CHIKV_MacB-x2352+C+401+1": 5962,
    "CHIKV_MacB-x1123_C_405_1_CHIKV_MacB-x0692+D+304+1": 5965,
    "CHIKV_MacB-x1123_D_306_1_CHIKV_MacB-x2352+C+401+1": 5963,
    "CHIKV_MacB-x1123_D_307_1_CHIKV_MacB-x2352+C+401+1": 5969,
    "CHIKV_MacB-x1125_A_401_1_CHIKV_MacB-x2352+C+401+1": 5972,
    "CHIKV_MacB-x1125_A_402_1_CHIKV_MacB-x1125+A+402+1": 5970,
    "CHIKV_MacB-x1125_C_401_1_CHIKV_MacB-x1125+C+401+1": 5971,
    "CHIKV_MacB-x1125_C_402_1_CHIKV_MacB-x2352+C+401+1": 5973,
    "CHIKV_MacB-x1131_A_401_1_CHIKV_MacB-x2352+C+401+1": 5977,
    "CHIKV_MacB-x1131_A_402_1_6w8q+A+201+1": 5974,
    "CHIKV_MacB-x1131_B_401_1_6w8q+A+201+1": 5975,
    "CHIKV_MacB-x1131_B_401_CHIKV_MacB-x0300+A+401+1": 75,
    "CHIKV_MacB-x1131_B_402_1_CHIKV_MacB-x2352+C+401+1": 5978,
    "CHIKV_MacB-x1131_B_402_CHIKV_MacB-x0300+A+401+1": 69,
    "CHIKV_MacB-x1131_C_308_1_CHIKV_MacB-x2352+C+401+1": 5979,
    "CHIKV_MacB-x1131_D_306_1_6w8q+A+201+1": 5976,
    "CHIKV_MacB-x1131_D_307_1_CHIKV_MacB-x2352+C+401+1": 5980,
    "CHIKV_MacB-x1132_A_305_1_CHIKV_MacB-x2352+C+401+1": 5981,
    "CHIKV_MacB-x1132_B_401_1_CHIKV_MacB-x2352+C+401+1": 5982,
    "CHIKV_MacB-x1132_B_401_CHIKV_MacB-x0300+A+401+1": 77,
    "CHIKV_MacB-x1132_B_403_1_CHIKV_MacB-x0543+B+308+1": 5986,
    "CHIKV_MacB-x1132_B_403_CHIKV_MacB-x1132+B+403+1": 88,
    "CHIKV_MacB-x1132_C_308_1_CHIKV_MacB-x2352+C+401+1": 5987,
    "CHIKV_MacB-x1132_C_308_CHIKV_MacB-x0300+A+401+1": 89,
    "CHIKV_MacB-x1132_D_304_1_CHIKV_MacB-x2352+C+401+1": 5983,
    "CHIKV_MacB-x1132_D_304_CHIKV_MacB-x0300+A+401+1": 76,
    "CHIKV_MacB-x1132_D_306_1_CHIKV_MacB-x2352+C+401+1": 5984,
    "CHIKV_MacB-x1132_D_307_1_CHIKV_MacB-x2352+C+401+1": 5985,
    "CHIKV_MacB-x1134_A_304_1_CHIKV_MacB-x1134+A+304+1": 5988,
    "CHIKV_MacB-x1134_A_305_1_CHIKV_MacB-x2352+C+401+1": 5991,
    "CHIKV_MacB-x1134_B_308_1_CHIKV_MacB-x1134+B+308+1": 5989,
    "CHIKV_MacB-x1134_B_309_1_CHIKV_MacB-x1134+B+309+1": 5990,
    "CHIKV_MacB-x1134_C_309_1_CHIKV_MacB-x2352+C+401+1": 5992,
    "CHIKV_MacB-x1134_D_401_1_CHIKV_MacB-x2352+C+401+1": 5993,
    "CHIKV_MacB-x1135_A_304_1_CHIKV_MacB-x2352+C+401+1": 5994,
    "CHIKV_MacB-x1135_A_305_1_CHIKV_MacB-x2352+C+401+1": 5995,
    "CHIKV_MacB-x1135_A_305_CHIKV_MacB-x0300+A+401+1": 37,
    "CHIKV_MacB-x1135_B_308_1_CHIKV_MacB-x2352+C+401+1": 5996,
    "CHIKV_MacB-x1135_B_308_CHIKV_MacB-x0300+A+401+1": 57,
    "CHIKV_MacB-x1135_B_309_1_CHIKV_MacB-x2352+C+401+1": 5997,
    "CHIKV_MacB-x1135_B_309_CHIKV_MacB-x0300+A+401+1": 58,
    "CHIKV_MacB-x1135_B_310_1_CHIKV_MacB-x2352+C+401+1": 6004,
    "CHIKV_MacB-x1135_B_310_CHIKV_MacB-x0300+A+401+1": 36,
    "CHIKV_MacB-x1135_B_593_1_CHIKV_MacB-x1134+B+308+1": 6003,
    "CHIKV_MacB-x1135_C_401_1_CHIKV_MacB-x2352+C+401+1": 5998,
    "CHIKV_MacB-x1135_C_403_1_CHIKV_MacB-x2352+C+401+1": 5999,
    "CHIKV_MacB-x1135_C_404_1_CHIKV_MacB-x2352+C+401+1": 6000,
    "CHIKV_MacB-x1135_D_401_1_CHIKV_MacB-x2352+C+401+1": 6005,
    "CHIKV_MacB-x1135_D_401_CHIKV_MacB-x0300+A+401+1": 51,
    "CHIKV_MacB-x1135_D_501_1_CHIKV_MacB-x2352+C+401+1": 6001,
    "CHIKV_MacB-x1135_D_501_CHIKV_MacB-x0300+A+401+1": 39,
    "CHIKV_MacB-x1135_D_502_1_CHIKV_MacB-x2352+C+401+1": 6002,
    "CHIKV_MacB-x1135_D_502_CHIKV_MacB-x0300+A+401+1": 40,
    "CHIKV_MacB-x1139_A_305_1_CHIKV_MacB-x2352+C+401+1": 5952,
    "CHIKV_MacB-x1139_C_306_1_CHIKV_MacB-x2352+C+401+1": 5953,
    "CHIKV_MacB-x1145_A_305_1_CHIKV_MacB-x2352+C+401+1": 6012,
    "CHIKV_MacB-x1145_B_308_1_CHIKV_MacB-x2352+C+401+1": 6006,
    "CHIKV_MacB-x1145_C_304_1_6w8q+A+201+1": 6009,
    "CHIKV_MacB-x1145_C_308_1_CHIKV_MacB-x0692+D+304+1": 6010,
    "CHIKV_MacB-x1145_C_309_1_CHIKV_MacB-x1123+B+310+1": 6011,
    "CHIKV_MacB-x1145_C_311_1_CHIKV_MacB-x2352+C+401+1": 6007,
    "CHIKV_MacB-x1145_D_306_1_CHIKV_MacB-x2352+C+401+1": 6008,
    "CHIKV_MacB-x1151_A_304_1_CHIKV_MacB-x2352+C+401+1": 6013,
    "CHIKV_MacB-x1151_A_304_CHIKV_MacB-x0300+A+401+1": 84,
    "CHIKV_MacB-x1151_A_305_1_CHIKV_MacB-x2352+C+401+1": 6014,
    "CHIKV_MacB-x1151_A_306_1_CHIKV_MacB-x2352+C+401+1": 6021,
    "CHIKV_MacB-x1151_B_308_1_CHIKV_MacB-x2352+C+401+1": 6020,
    "CHIKV_MacB-x1151_B_309_1_CHIKV_MacB-x2352+C+401+1": 6015,
    "CHIKV_MacB-x1151_C_401_1_CHIKV_MacB-x2352+C+401+1": 6016,
    "CHIKV_MacB-x1151_C_401_CHIKV_MacB-x0300+A+401+1": 86,
    "CHIKV_MacB-x1151_C_402_1_CHIKV_MacB-x2352+C+401+1": 6017,
    "CHIKV_MacB-x1151_C_403_1_CHIKV_MacB-x2352+C+401+1": 6023,
    "CHIKV_MacB-x1151_C_403_CHIKV_MacB-x0300+A+401+1": 82,
    "CHIKV_MacB-x1151_D_401_1_CHIKV_MacB-x2352+C+401+1": 6022,
    "CHIKV_MacB-x1151_D_402_1_CHIKV_MacB-x2352+C+401+1": 6018,
    "CHIKV_MacB-x1151_D_403_1_CHIKV_MacB-x2352+C+401+1": 6019,
    "CHIKV_MacB-x1175_A_401_1_CHIKV_MacB-x2352+C+401+1": 6024,
    "CHIKV_MacB-x1175_B_304_1_CHIKV_MacB-x2352+C+401+1": 6025,
    "CHIKV_MacB-x1182_A_305_1_CHIKV_MacB-x2352+C+401+1": 6026,
    "CHIKV_MacB-x1182_A_305_CHIKV_MacB-x0300+A+401+1": 14,
    "CHIKV_MacB-x1182_B_401_1_CHIKV_MacB-x2352+C+401+1": 6027,
    "CHIKV_MacB-x1182_D_306_1_CHIKV_MacB-x2352+C+401+1": 6028,
    "CHIKV_MacB-x1183_B_304_1_CHIKV_MacB-x2352+C+401+1": 6030,
    "CHIKV_MacB-x1183_B_401_1_CHIKV_MacB-x2352+C+401+1": 6029,
    "CHIKV_MacB-x1185_A_305_1_CHIKV_MacB-x2352+C+401+1": 6031,
    "CHIKV_MacB-x1185_B_401_1_CHIKV_MacB-x2352+C+401+1": 6032,
    "CHIKV_MacB-x1185_C_304_1_CHIKV_MacB-x2352+C+401+1": 6033,
    "CHIKV_MacB-x1185_D_401_1_CHIKV_MacB-x2352+C+401+1": 6034,
    "CHIKV_MacB-x1188_C_401_1_CHIKV_MacB-x2352+C+401+1": 6035,
    "CHIKV_MacB-x1194_A_305_1_CHIKV_MacB-x2352+C+401+1": 6036,
    "CHIKV_MacB-x1194_D_401_1_CHIKV_MacB-x2352+C+401+1": 6037,
    "CHIKV_MacB-x1198_A_401_1_CHIKV_MacB-x2352+C+401+1": 5954,
    "CHIKV_MacB-x1198_B_304_1_CHIKV_MacB-x2352+C+401+1": 5955,
    "CHIKV_MacB-x1198_C_308_1_CHIKV_MacB-x2352+C+401+1": 5956,
    "CHIKV_MacB-x1198_D_401_1_CHIKV_MacB-x2352+C+401+1": 5957,
    "CHIKV_MacB-x1203_A_304_1_CHIKV_MacB-x2352+C+401+1": 6040,
    "CHIKV_MacB-x1203_A_401_1_CHIKV_MacB-x2352+C+401+1": 6041,
    "CHIKV_MacB-x1203_A_401_CHIKV_MacB-x0300+A+401+1": 33,
    "CHIKV_MacB-x1203_A_501_1_CHIKV_MacB-x2352+C+401+1": 6042,
    "CHIKV_MacB-x1203_A_501_CHIKV_MacB-x0300+A+401+1": 10,
    "CHIKV_MacB-x1203_B_308_1_CHIKV_MacB-x2352+C+401+1": 6043,
    "CHIKV_MacB-x1203_B_308_CHIKV_MacB-x0300+A+401+1": 12,
    "CHIKV_MacB-x1203_B_309_1_CHIKV_MacB-x2352+C+401+1": 6045,
    "CHIKV_MacB-x1203_B_309_CHIKV_MacB-x0300+A+401+1": 30,
    "CHIKV_MacB-x1203_C_309_1_CHIKV_MacB-x2352+C+401+1": 6046,
    "CHIKV_MacB-x1203_C_310_1_CHIKV_MacB-x2352+C+401+1": 6044,
    "CHIKV_MacB-x1203_D_306_1_CHIKV_MacB-x2352+C+401+1": 6047,
    "CHIKV_MacB-x1213_A_401_1_CHIKV_MacB-x2352+C+401+1": 6048,
    "CHIKV_MacB-x1213_A_401_CHIKV_MacB-x0300+A+401+1": 90,
    "CHIKV_MacB-x1213_A_402_1_CHIKV_MacB-x2352+C+401+1": 6049,
    "CHIKV_MacB-x1213_B_308_1_CHIKV_MacB-x2352+C+401+1": 6050,
    "CHIKV_MacB-x1213_C_401_1_CHIKV_MacB-x2352+C+401+1": 6051,
    "CHIKV_MacB-x1213_D_304_1_CHIKV_MacB-x2352+C+401+1": 6052,
    "CHIKV_MacB-x1258_A_305_1_CHIKV_MacB-x2352+C+401+1": 6053,
    "CHIKV_MacB-x1258_B_304_1_CHIKV_MacB-x2352+C+401+1": 6054,
    "CHIKV_MacB-x1258_C_308_1_CHIKV_MacB-x2352+C+401+1": 6055,
    "CHIKV_MacB-x1258_D_401_1_CHIKV_MacB-x2352+C+401+1": 6056,
    "CHIKV_MacB-x1260_A_304_1_CHIKV_MacB-x2352+C+401+1": 6059,
    "CHIKV_MacB-x1260_A_304_CHIKV_MacB-x0300+A+401+1": 21,
    "CHIKV_MacB-x1260_A_305_1_CHIKV_MacB-x2352+C+401+1": 6057,
    "CHIKV_MacB-x1260_B_401_1_CHIKV_MacB-x2352+C+401+1": 6058,
    "CHIKV_MacB-x1268_B_304_1_CHIKV_MacB-x2352+C+401+1": 6060,
    "CHIKV_MacB-x1268_C_401_1_CHIKV_MacB-x2352+C+401+1": 6061,
    "CHIKV_MacB-x1268_D_306_1_CHIKV_MacB-x2352+C+401+1": 6062,
    "CHIKV_MacB-x1271_A_401_1_CHIKV_MacB-x2352+C+401+1": 6063,
    "CHIKV_MacB-x1271_A_401_CHIKV_MacB-x0300+A+401+1": 38,
    "CHIKV_MacB-x1271_B_308_1_CHIKV_MacB-x2352+C+401+1": 6064,
    "CHIKV_MacB-x1271_C_401_1_CHIKV_MacB-x2352+C+401+1": 6065,
    "CHIKV_MacB-x1271_C_402_1_CHIKV_MacB-x2352+C+401+1": 6066,
    "CHIKV_MacB-x1271_D_304_1_CHIKV_MacB-x2352+C+401+1": 6068,
    "CHIKV_MacB-x1271_D_306_1_CHIKV_MacB-x2352+C+401+1": 6067,
    "CHIKV_MacB-x1271_D_306_CHIKV_MacB-x0300+A+401+1": 35,
    "CHIKV_MacB-x1278_D_401_1_CHIKV_MacB-x2352+C+401+1": 6069,
    "CHIKV_MacB-x1278_D_401_CHIKV_MacB-x0300+A+401+1": 67,
    "CHIKV_MacB-x1288_A_401_1_CHIKV_MacB-x2352+C+401+1": 6070,
    "CHIKV_MacB-x1288_A_402_1_CHIKV_MacB-x2352+C+401+1": 6071,
    "CHIKV_MacB-x1288_A_403_1_CHIKV_MacB-x2352+C+401+1": 6075,
    "CHIKV_MacB-x1288_B_304_1_CHIKV_MacB-x2352+C+401+1": 6072,
    "CHIKV_MacB-x1288_D_401_1_CHIKV_MacB-x2352+C+401+1": 6073,
    "CHIKV_MacB-x1288_D_402_1_CHIKV_MacB-x2352+C+401+1": 6074,
    "CHIKV_MacB-x1298_B_401_1_CHIKV_MacB-x2352+C+401+1": 6038,
    "CHIKV_MacB-x1298_B_401_CHIKV_MacB-x0300+A+401+1": 41,
    "CHIKV_MacB-x1298_D_401_1_CHIKV_MacB-x2352+C+401+1": 6039,
    "CHIKV_MacB-x1299_B_304_1_CHIKV_MacB-x1123+B+310+1": 6077,
    "CHIKV_MacB-x1299_B_304_CHIKV_MacB-x1123+B+310+1": 81,
    "CHIKV_MacB-x1299_C_401_1_CHIKV_MacB-x2352+C+401+1": 6078,
    "CHIKV_MacB-x1299_D_401_1_CHIKV_MacB-x2352+C+401+1": 6079,
    "CHIKV_MacB-x1305_A_304_1_CHIKV_MacB-x2352+C+401+1": 6080,
    "CHIKV_MacB-x1305_A_304_CHIKV_MacB-x0300+A+401+1": 48,
    "CHIKV_MacB-x1305_A_305_1_CHIKV_MacB-x2352+C+401+1": 6081,
    "CHIKV_MacB-x1305_B_308_1_CHIKV_MacB-x2352+C+401+1": 6082,
    "CHIKV_MacB-x1305_C_401_1_CHIKV_MacB-x2352+C+401+1": 6083,
    "CHIKV_MacB-x1305_D_306_1_CHIKV_MacB-x2352+C+401+1": 6084,
    "CHIKV_MacB-x1315_C_304_1_CHIKV_MacB-x2352+C+401+1": 6085,
    "CHIKV_MacB-x1315_D_306_1_CHIKV_MacB-x2352+C+401+1": 6086,
    "CHIKV_MacB-x1316_A_401_1_CHIKV_MacB-x2352+C+401+1": 6088,
    "CHIKV_MacB-x1316_C_304_1_CHIKV_MacB-x1316+C+304+1": 6087,
    "CHIKV_MacB-x1316_D_306_1_CHIKV_MacB-x2352+C+401+1": 6089,
    "CHIKV_MacB-x1320_D_304_1_CHIKV_MacB-x2352+C+401+1": 5798,
    "CHIKV_MacB-x1320_D_304_CHIKV_MacB-x0300+A+401+1": 61,
    "CHIKV_MacB-x1331_A_304_1_CHIKV_MacB-x1134+B+308+1": 6090,
    "CHIKV_MacB-x1331_A_305_1_CHIKV_MacB-x2352+C+401+1": 6091,
    "CHIKV_MacB-x1331_B_308_1_CHIKV_MacB-x2352+C+401+1": 6092,
    "CHIKV_MacB-x1331_C_401_1_CHIKV_MacB-x2352+C+401+1": 6093,
    "CHIKV_MacB-x1331_D_306_1_CHIKV_MacB-x2352+C+401+1": 6094,
    "CHIKV_MacB-x1335_B_401_1_CHIKV_MacB-x2352+C+401+1": 6095,
    "CHIKV_MacB-x1335_C_304_1_CHIKV_MacB-x2352+C+401+1": 6096,
    "CHIKV_MacB-x1335_D_401_1_CHIKV_MacB-x2352+C+401+1": 6097,
    "CHIKV_MacB-x1337_C_401_1_CHIKV_MacB-x2352+C+401+1": 6099,
    "CHIKV_MacB-x1337_D_304_1_CHIKV_MacB-x0692+D+304+1": 6098,
    "CHIKV_MacB-x1337_D_401_1_CHIKV_MacB-x2352+C+401+1": 6100,
    "CHIKV_MacB-x1338_B_304_1_CHIKV_MacB-x2352+C+401+1": 6101,
    "CHIKV_MacB-x1338_B_304_CHIKV_MacB-x0300+A+401+1": 18,
    "CHIKV_MacB-x1338_C_308_1_CHIKV_MacB-x2352+C+401+1": 6102,
    "CHIKV_MacB-x1338_C_309_1_CHIKV_MacB-x2352+C+401+1": 6103,
    "CHIKV_MacB-x1353_D_304_1_CHIKV_MacB-x2352+C+401+1": 6106,
    "CHIKV_MacB-x1353_D_401_1_CHIKV_MacB-x2352+C+401+1": 6104,
    "CHIKV_MacB-x1353_D_501_1_CHIKV_MacB-x2352+C+401+1": 6105,
    "CHIKV_MacB-x1365_A_304_1_CHIKV_MacB-x2352+C+401+1": 6107,
    "CHIKV_MacB-x1365_A_401_1_CHIKV_MacB-x2352+C+401+1": 6108,
    "CHIKV_MacB-x1365_A_401_CHIKV_MacB-x0300+A+401+1": 42,
    "CHIKV_MacB-x1365_C_308_1_CHIKV_MacB-x2352+C+401+1": 6109,
    "CHIKV_MacB-x1365_C_308_CHIKV_MacB-x0300+A+401+1": 66,
    "CHIKV_MacB-x1371_C_401_1_CHIKV_MacB-x2352+C+401+1": 6112,
    "CHIKV_MacB-x1371_C_401_CHIKV_MacB-x0300+A+401+1": 19,
    "CHIKV_MacB-x1371_C_402_1_CHIKV_MacB-x2352+C+401+1": 6110,
    "CHIKV_MacB-x1371_C_403_1_CHIKV_MacB-x2352+C+401+1": 6111,
    "CHIKV_MacB-x1372_A_304_1_CHIKV_MacB-x2352+C+401+1": 6113,
    "CHIKV_MacB-x1372_A_304_CHIKV_MacB-x0300+A+401+1": 4,
    "CHIKV_MacB-x1372_B_308_1_CHIKV_MacB-x2352+C+401+1": 6114,
    "CHIKV_MacB-x1372_C_401_1_CHIKV_MacB-x2352+C+401+1": 6115,
    "CHIKV_MacB-x1372_D_306_1_CHIKV_MacB-x2352+C+401+1": 6116,
    "CHIKV_MacB-x1379_A_305_1_CHIKV_MacB-x2352+C+401+1": 6117,
    "CHIKV_MacB-x1379_A_306_1_CHIKV_MacB-x2352+C+401+1": 6118,
    "CHIKV_MacB-x1379_B_308_1_CHIKV_MacB-x2352+C+401+1": 6122,
    "CHIKV_MacB-x1379_C_401_1_CHIKV_MacB-x2352+C+401+1": 6123,
    "CHIKV_MacB-x1379_C_501_1_CHIKV_MacB-x2352+C+401+1": 6119,
    "CHIKV_MacB-x1379_D_304_1_CHIKV_MacB-x2352+C+401+1": 6120,
    "CHIKV_MacB-x1379_D_401_1_CHIKV_MacB-x2352+C+401+1": 6121,
    "CHIKV_MacB-x1397_A_304_1_CHIKV_MacB-x2352+C+401+1": 6124,
    "CHIKV_MacB-x1397_D_401_1_CHIKV_MacB-x2352+C+401+1": 6125,
    "CHIKV_MacB-x1402_A_401_1_CHIKV_MacB-x2352+C+401+1": 6126,
    "CHIKV_MacB-x1402_A_501_1_CHIKV_MacB-x2352+C+401+1": 6127,
    "CHIKV_MacB-x1402_A_501_CHIKV_MacB-x0300+A+401+1": 44,
    "CHIKV_MacB-x1402_B_304_1_CHIKV_MacB-x2352+C+401+1": 6128,
    "CHIKV_MacB-x1402_B_304_CHIKV_MacB-x0300+A+401+1": 64,
    "CHIKV_MacB-x1402_B_308_1_CHIKV_MacB-x2352+C+401+1": 6129,
    "CHIKV_MacB-x1402_B_308_CHIKV_MacB-x0300+A+401+1": 62,
    "CHIKV_MacB-x1402_C_308_1_CHIKV_MacB-x2352+C+401+1": 6131,
    "CHIKV_MacB-x1402_C_308_CHIKV_MacB-x0300+A+401+1": 46,
    "CHIKV_MacB-x1402_C_309_1_CHIKV_MacB-x2352+C+401+1": 6130,
    "CHIKV_MacB-x1402_D_401_1_CHIKV_MacB-x2352+C+401+1": 6132,
    "CHIKV_MacB-x1406_A_304_1_CHIKV_MacB-x2352+C+401+1": 6133,
    "CHIKV_MacB-x1406_B_401_1_CHIKV_MacB-x2352+C+401+1": 6134,
    "CHIKV_MacB-x1410_A_401_1_CHIKV_MacB-x2352+C+401+1": 6135,
    "CHIKV_MacB-x1410_B_304_1_CHIKV_MacB-x2352+C+401+1": 6136,
    "CHIKV_MacB-x1410_B_304_CHIKV_MacB-x0300+A+401+1": 11,
    "CHIKV_MacB-x1410_C_308_1_CHIKV_MacB-x2352+C+401+1": 6137,
    "CHIKV_MacB-x1410_C_308_CHIKV_MacB-x0300+A+401+1": 32,
    "CHIKV_MacB-x1410_D_401_1_CHIKV_MacB-x2352+C+401+1": 6138,
    "CHIKV_MacB-x1410_D_401_CHIKV_MacB-x0300+A+401+1": 31,
    "CHIKV_MacB-x1417_A_304_1_CHIKV_MacB-x2352+C+401+1": 6139,
    "CHIKV_MacB-x1417_B_308_1_CHIKV_MacB-x2352+C+401+1": 6140,
    "CHIKV_MacB-x1417_C_308_1_CHIKV_MacB-x2352+C+401+1": 6141,
    "CHIKV_MacB-x1417_D_401_1_CHIKV_MacB-x2352+C+401+1": 6142,
    "CHIKV_MacB-x1421_A_401_1_CHIKV_MacB-x2352+C+401+1": 6143,
    "CHIKV_MacB-x1421_C_304_1_CHIKV_MacB-x2352+C+401+1": 6144,
    "CHIKV_MacB-x1423_C_304_1_CHIKV_MacB-x2352+C+401+1": 6076,
    "CHIKV_MacB-x1431_B_304_1_CHIKV_MacB-x2352+C+401+1": 6145,
    "CHIKV_MacB-x1431_B_308_1_CHIKV_MacB-x2352+C+401+1": 6146,
    "CHIKV_MacB-x1431_D_401_1_CHIKV_MacB-x2352+C+401+1": 6147,
    "CHIKV_MacB-x1431_D_402_1_CHIKV_MacB-x2352+C+401+1": 6148,
    "CHIKV_MacB-x1435_C_401_1_CHIKV_MacB-x2352+C+401+1": 6149,
    "CHIKV_MacB-x1440_A_305_1_CHIKV_MacB-x2352+C+401+1": 6150,
    "CHIKV_MacB-x1440_A_306_1_CHIKV_MacB-x2352+C+401+1": 6151,
    "CHIKV_MacB-x1440_B_401_1_CHIKV_MacB-x2352+C+401+1": 6152,
    "CHIKV_MacB-x1440_C_401_1_CHIKV_MacB-x2352+C+401+1": 6153,
    "CHIKV_MacB-x1440_C_501_1_CHIKV_MacB-x2352+C+401+1": 6155,
    "CHIKV_MacB-x1440_D_304_1_CHIKV_MacB-x2352+C+401+1": 6154,
    "CHIKV_MacB-x1440_D_401_1_CHIKV_MacB-x2352+C+401+1": 6156,
    "CHIKV_MacB-x1440_D_401_CHIKV_MacB-x0300+A+401+1": 85,
    "CHIKV_MacB-x1444_A_401_1_CHIKV_MacB-x2352+C+401+1": 6157,
    "CHIKV_MacB-x1444_A_401_CHIKV_MacB-x0300+A+401+1": 45,
    "CHIKV_MacB-x1444_B_304_1_CHIKV_MacB-x2352+C+401+1": 6158,
    "CHIKV_MacB-x1444_D_306_1_CHIKV_MacB-x2352+C+401+1": 6159,
    "CHIKV_MacB-x1445_A_305_1_CHIKV_MacB-x2352+C+401+1": 6160,
    "CHIKV_MacB-x1445_A_306_1_CHIKV_MacB-x2352+C+401+1": 6161,
    "CHIKV_MacB-x1445_B_401_1_CHIKV_MacB-x2352+C+401+1": 6162,
    "CHIKV_MacB-x1445_B_401_CHIKV_MacB-x0300+A+401+1": 73,
    "CHIKV_MacB-x1445_B_402_1_CHIKV_MacB-x2352+C+401+1": 6163,
    "CHIKV_MacB-x1445_B_402_CHIKV_MacB-x0300+A+401+1": 70,
    "CHIKV_MacB-x1445_C_401_1_CHIKV_MacB-x2352+C+401+1": 6165,
    "CHIKV_MacB-x1445_D_304_1_CHIKV_MacB-x2352+C+401+1": 6164,
    "CHIKV_MacB-x1454_A_401_1_CHIKV_MacB-x2352+C+401+1": 6166,
    "CHIKV_MacB-x1454_B_401_1_CHIKV_MacB-x2352+C+401+1": 6167,
    "CHIKV_MacB-x1454_C_308_1_CHIKV_MacB-x2352+C+401+1": 6168,
    "CHIKV_MacB-x1454_D_304_1_CHIKV_MacB-x2352+C+401+1": 6169,
    "CHIKV_MacB-x2129_A_304_1_CHIKV_MacB-x2352+C+401+1": 6170,
    "CHIKV_MacB-x2129_B_308_1_CHIKV_MacB-x2352+C+401+1": 6171,
    "CHIKV_MacB-x2129_C_308_1_CHIKV_MacB-x2352+C+401+1": 6172,
    "CHIKV_MacB-x2129_D_306_1_CHIKV_MacB-x2352+C+401+1": 6173,
    "CHIKV_MacB-x2182_A_304_1_CHIKV_MacB-x2352+C+401+1": 6174,
    "CHIKV_MacB-x2182_B_204_1_CHIKV_MacB-x2352+C+401+1": 6175,
    "CHIKV_MacB-x2182_C_308_1_CHIKV_MacB-x2352+C+401+1": 6176,
    "CHIKV_MacB-x2182_D_206_1_CHIKV_MacB-x2352+C+401+1": 6177,
    "CHIKV_MacB-x2183_A_304_1_CHIKV_MacB-x2352+C+401+1": 6178,
    "CHIKV_MacB-x2183_B_205_1_CHIKV_MacB-x2352+C+401+1": 6179,
    "CHIKV_MacB-x2183_C_308_1_CHIKV_MacB-x2352+C+401+1": 6180,
    "CHIKV_MacB-x2183_D_206_1_CHIKV_MacB-x2352+C+401+1": 6181,
    "CHIKV_MacB-x2343_C_401_2_CHIKV_MacB-x2346+C+401+1": 5797,
    "CHIKV_MacB-x2344_B_401_1_CHIKV_MacB-x2352+C+401+1": 6182,
    "CHIKV_MacB-x2346_C_401_1_CHIKV_MacB-x2346+C+401+1": 6183,
    "CHIKV_MacB-x2350_A_501_1_CHIKV_MacB-x2352+C+401+1": 6184,
    "CHIKV_MacB-x2350_B_308_1_CHIKV_MacB-x2352+C+401+1": 6185,
    "CHIKV_MacB-x2350_D_306_1_CHIKV_MacB-x2352+C+401+1": 6186,
    "CHIKV_MacB-x2352_C_401_1_CHIKV_MacB-x2352+C+401+1": 6187,
    "CHIKV_MacB-x2364_A_304_1_CHIKV_MacB-x2352+C+401+1": 6188,
    "CHIKV_MacB-x2364_B_401_1_CHIKV_MacB-x2352+C+401+1": 6189,
    "CHIKV_MacB-x2364_D_401_1_CHIKV_MacB-x2352+C+401+1": 6190,
    "CHIKV_MacB-x2370_B_401_2_CHIKV_MacB-x2352+C+401+1": 6191,
    "CHIKV_MacB-x2431_B_308_2_CHIKV_MacB-x2352+C+401+1": 6192,
}

if __name__ == "__main__":
    unittest.main()
