import unittest
import os, sys
sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals


# TODO: test composite that changes number of entries

class composite(unittest.TestCase):
    def assertArraysAlmostEqual(self,a,b):
        self.assertEqual(len(a),len(b))
        for (i1,i2) in zip(a,b):
            self.assertAlmostEqual(i1,i2)

    def test_unroll(self):

        min1 = minerals.Murakami_etal_2012.fe_periclase()
        min_hs = burnman.minerals.Murakami_etal_2012.fe_periclase_HS()
        min2 = minerals.SLB_2005.periclase()

        print min1
        c = burnman.composite( [(min1,1.0)] )
        c.set_method("slb2")
        c.set_state(5e9,300)
        (f,m) = c.unroll()
        mins=",".join([a.to_string() for a in m])
        self.assertEqual(f,[1.0])
        self.assertEqual(mins,min_hs.to_string())
        c = burnman.composite( [(min1,0.4),(min2,0.6)] )
        c.set_method("slb2")
        c.set_state(5e9,300)
        (f,m) = c.unroll()
        self.assertEqual(f,[0.4,0.6])

        c1 = burnman.composite( [(min1,1.0)] )
        c2 = burnman.composite( [(min2,1.0)] )
        c = burnman.composite( [(min1,0.1),(c1,0.4),(c2,0.5)] )
        (f,m) = c.unroll()
        mins=",".join([a.to_string() for a in m])
        self.assertEqual(f,[0.1,0.4,0.5])
        self.assertEqual(mins,min_hs.to_string()+','+min_hs.to_string()+','+min2.to_string())

        c1 = burnman.composite( [ (min1,0.1), (min2,0.9) ] )
        c2 = burnman.composite( [(min2,1.0)] )
        c = burnman.composite( [(min1,0.3),(c1,0.1),(c2,0.6)] )
        (f,m) = c.unroll()
        mins=",".join([a.to_string() for a in m])
        self.assertArraysAlmostEqual(f,[0.3,0.01,0.09,0.6])
        self.assertEqual(mins,",".join([min_hs.to_string(),min_hs.to_string(),min2.to_string(),min2.to_string()]))
        
        #c.set_method("slb2")
        #c.set_state(5e9,300)

    def test_changevalues(self):
        class mycomposite(burnman.material):
            def unroll(self):
                fractions = [0.3, 0.7]
                mins = [minerals.Murakami_etal_2012.fe_periclase(), minerals.SLB_2005.periclase()]
                if (self.temperature>500):
                    fractions = [0.1, 0.9]
                return (fractions, mins)
        
        c = mycomposite()
        c.set_state(5e9,300)
        (f,m) = c.unroll()
        mins=",".join([a.to_string() for a in m])
        self.assertArraysAlmostEqual(f,[0.3,0.7])
        self.assertEqual(mins,",".join([minerals.Murakami_etal_2012.fe_periclase().to_string(),minerals.SLB_2005.periclase().to_string()]))
        c.set_state(5e9,3000)
        (f,m) = c.unroll()
        self.assertArraysAlmostEqual(f,[0.1,0.9])
        self.assertEqual(mins,",".join([minerals.Murakami_etal_2012.fe_periclase().to_string(),minerals.SLB_2005.periclase().to_string()]))

    def test_number(self):
        class mycomposite(burnman.material):
            def unroll(self):
                if (self.temperature>500):
                    return ([1.0],[minerals.Murakami_etal_2012.fe_periclase()])
                fractions = [0.3, 0.7]
                mins = [minerals.Murakami_etal_2012.fe_periclase(), minerals.SLB_2005.periclase()]
                return (fractions, mins)
        
        c = mycomposite()
        c.set_state(5e9,1000)
        (f,m) = c.unroll()
        mins=",".join([a.to_string() for a in m])
        self.assertArraysAlmostEqual(f,[1.0])
        self.assertEqual(mins,minerals.Murakami_etal_2012.fe_periclase().to_string())
        c.set_state(5e9,300)
        (f,m) = c.unroll()
        mins=",".join([a.to_string() for a in m])
        self.assertArraysAlmostEqual(f,[0.3,0.7])
        self.assertEqual(mins,",".join([minerals.Murakami_etal_2012.fe_periclase().to_string(),minerals.SLB_2005.periclase().to_string()]))

    def test_nest(self):
        min1 = minerals.Murakami_etal_2012.fe_periclase_LS()
        min2 = minerals.SLB_2005.periclase()
        ca = burnman.composite( [(min1,1.0)] )
        c = burnman.composite( [(ca,0.4),(min2, 0.6)] )
        c.set_method("slb3")
        c.set_state(5e9,1000)
        (f,m) = c.unroll()
        mins=",".join([a.to_string() for a in m])
        self.assertArraysAlmostEqual(f,[0.4,0.6])
        self.assertEqual(mins,",".join([min1.to_string(),min2.to_string()]))

if __name__ == '__main__':
    unittest.main()
