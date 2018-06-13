# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 09:11:42 2018

@author: yiyuezhuo
"""

import unittest
from core import OneSampleTTestPower,TwoSampleTTestPower

class TestOneSampleTTestPower(unittest.TestCase):
    
    def test_reject_region(self):
        '''
        The result in G*power is not expoesd to user.
        For essential, the result is due only to alpha and statistic
        distribution specification in H0.
        '''
        test = OneSampleTTestPower(effect_size = 0.5,alpha=0.05,power=0.95)
        
        test.N = 10
        test.reject_region = test.get_reject_region()
        
        self.assertEqual(test.reject_region['type'],'le')
        self.assertAlmostEqual(test.reject_region['value'],1.8331129326536335)
        '''
              ___
             /   \|
            /     \
           /      |\
          /       | \
         /        |  \
        /         |5% \
              0  1.88 
        '''
        
        test.N = 20
        test.reject_region = test.get_reject_region()
        
        self.assertEqual(test.reject_region['type'],'le')
        self.assertAlmostEqual(test.reject_region['value'],1.729132811521367)


    def test_get_power(self):
        '''
        G*Power---Means: Difference from constant (one sample case)---Pose hoc
        Effect size d = 0.5
        \alpha error prob = 0.05
        Total sample size = 100
        '''
        test = OneSampleTTestPower(effect_size = 0.5,alpha=0.05,power=0.95)
        test.N = 100
        test.reject_region = test.get_reject_region()
        power = test.get_power()
        #self.assertAlmostEqual(power, 0.999550905486557) # python result
        self.assertAlmostEqual(power, 0.9995509) # G*Power result
        
    def test_get_N(self):
        '''
        G*Power---Means: Difference from constant (one sample case)---A priori
        Effect size d = 0.5
        \alpha error prob = 0.05
        Power(1 - \beta error prob) = 0.95
        '''
        test = OneSampleTTestPower(effect_size = 0.5,alpha=0.05,power=0.95)
        N = test.get_N()
        self.assertEqual(N, 45)

class TestTwoSampleTTestPower(unittest.TestCase):
    
    def test_get_N(self):
        test = TwoSampleTTestPower(effect_size = 0.5,alpha=0.05,power=0.95)
        N = test.get_N()
        self.assertEqual(N,176)
        
    def test_get_power(self):
        
        test = TwoSampleTTestPower(effect_size = 0.5,alpha=0.05,power=0.95)
        test.N = 30
        test.reject_region = test.get_reject_region()
        self.assertAlmostEqual(test.get_power(),0.379219288756931)
        
        test = TwoSampleTTestPower(effect_size = 0.5,alpha=0.05,power=0.95)
        test.N = 175
        test.reject_region = test.get_reject_region()
        self.assertAlmostEqual(test.get_power(),0.9494984955069458)

        test = TwoSampleTTestPower(effect_size = 0.5,alpha=0.05,power=0.95)
        test.N = 176
        test.reject_region = test.get_reject_region()
        self.assertAlmostEqual(test.get_power(),0.9514329065971822)

if __name__ == '__main__':
    unittest.main()
