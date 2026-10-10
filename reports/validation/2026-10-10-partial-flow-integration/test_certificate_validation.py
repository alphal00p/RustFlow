"""Synthetic schema/binding rejection tests, not tests of the physical theorem."""
import copy
import unittest
from certificate_validation import *

CUTS=[[],[0],[1],[0,1]]
def fixture():
    definition={'loops':2,'edges':[{'mass_squared':'0'}for _ in range(3)],'targets':[{},{}]}
    records=[]
    for cut in [0,1]:
        shifted=[i for i in range(3)if i!=cut]
        audit={'origin_identity':ORIGIN,'endpoint_structure':{'proof_version':ENDPOINT,'virtual_loops':1,
            'denominator_rows':[[i,['1'],'1']for i in shifted],'block_assignments':[],'block_representatives':[],
            'block_a_lower_bounds':[],'block_q_upper_bounds':['1']},
            'labels':[{'indices':[1,1,1,0,0,0,0],'classification':'uniform-null-cone-endpoint',
            'high_dimension_witness':'13/2','inequalities':[{'degree':{'twice_dimension_coefficient':1,'constant':-3}}]}]}
        records.append({'cut_slots':[cut],'certificate':{'construction':CONSTRUCTION,'proof_version':PROOF,
            'scope':'original eta-independent physical singleton targets only','native_closure':False,'transport_performed':False,
            'input_identity':'synthetic input','family_signature':[str(cut)],'proof_homotopy_identity':'synthetic proof',
            'normalization':'synthetic zero normalization','proof_homotopy_shifted_slots':shifted,'label_audit':audit,'nonzero_conditions':['epsilon+1']}})
    rows=[]
    for item in records:
        c=item['certificate'];rows.append({'cut_slots':item['cut_slots'],'report':{'construction':CONSTRUCTION,'basis_size':0,
            'native_storage_capacity':None,'source_options':None,'shifted_slots':[],'massless_endpoint':c['label_audit'],
            'nonzero_conditions':c['nonzero_conditions'],'contour_admission':PROOF,'physical_arity':7,'values':['0','(0+0i)']}})
    options={'policy':'synthetic'}
    rows.append({'cut_slots':[0,1],'report':{'construction':'weighted_amf','massless_endpoint':{'test':True},'source_options':options}})
    sample={'physical_zero_certificates':copy.deepcopy(records),'occupied_reports':rows,
            'contributions':[{'cut_slots':c,'values':['0','0']}for c in CUTS], 'vacuum_zero_certificates':[{'synthetic':True}],'empty_support':[]}
    return definition,records,sample,options

class CertificateValidation(unittest.TestCase):
    def test_well_formed_distinct_constructions(self):
        d,c,r,o=fixture();validate_certificates(c,d,CUTS);validate_sample(r,c,d,CUTS,o)
    def test_reject_bad_certificate_contracts(self):
        changes=[lambda c:c.pop(),lambda c:c.reverse(),lambda c:c[0]['certificate'].update(native_closure=True),
                 lambda c:c[0]['certificate'].update(transport_performed=True),lambda c:c[0]['certificate'].update(proof_version='other'),
                 lambda c:c[0]['certificate'].update(proof_homotopy_shifted_slots=[2]),
                 lambda c:c[0]['certificate']['label_audit']['endpoint_structure'].update(virtual_loops=0),
                 lambda c:c[0]['certificate']['label_audit']['endpoint_structure'].update(block_assignments=[[1,0]]),
                 lambda c:c[0]['certificate']['label_audit']['labels'][0].update(high_dimension_witness='6'),
                 lambda c:c[0]['certificate']['label_audit']['labels'][0].update(classification='independent-rank-one-block-endpoint')]
        for change in changes:
            with self.subTest(change=changes.index(change)):
                d,c,_,_=fixture();change(c)
                with self.assertRaises((ValueError,KeyError)):validate_certificates(c,d,CUTS)
    def test_reject_sample_binding_or_construction_changes(self):
        changes=[lambda r:r['physical_zero_certificates'][0]['certificate'].update(input_identity='different'),
            lambda r:r['occupied_reports'][0]['report'].update(construction='weighted_amf'),
            lambda r:r['occupied_reports'][0]['report'].update(basis_size=1),
            lambda r:r['occupied_reports'][0]['report'].update(native_storage_capacity=7),
            lambda r:r['occupied_reports'][0]['report'].update(source_options={'policy':'synthetic'}),
            lambda r:r['occupied_reports'][0]['report'].update(shifted_slots=[1]),
            lambda r:r['occupied_reports'][0]['report'].update(nonzero_conditions=[]),
            lambda r:r['occupied_reports'][0]['report'].update(values=['(1e-100+0i)','0']),
            lambda r:r['contributions'][1].update(values=['(1e-100+0i)','0']),
            lambda r:r['occupied_reports'][2]['report'].update(construction=CONSTRUCTION),
            lambda r:r['occupied_reports'][2]['report'].update(source_options={}),
            lambda r:r.update(vacuum_zero_certificates=[]),lambda r:r.update(empty_support=[[0]])]
        for change in changes:
            with self.subTest(change=changes.index(change)):
                d,c,r,o=fixture();change(r)
                with self.assertRaises((ValueError,KeyError)):validate_sample(r,c,d,CUTS,o)
if __name__=='__main__':unittest.main()
