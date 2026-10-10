import copy,unittest
from history_coverage import verify_history_coverage

def fixture(limit):
    history=[[0,1],[1,0],[2,0]];count=len(history)if limit is None else min(limit,len(history))
    prior={'status':'provisional','round':3,'active_target_closure':True,'storage_capacity':2,
        'max_history_candidate_maps':limit,'program_bytes':100,'program_blake3':'prior',
        'active_state':{'historical_requests':history,'submitted_requests':[[1,0]],
            'next_needed':[[0,1]],'history_is_coverage_certificate':False}}
    closed=copy.deepcopy(prior);closed.update(status='closed',program_bytes=104,program_blake3='closed')
    closed['active_state']['historical_candidate_collection']={
        'max_history_candidate_maps':limit,'selection':'ascending stored-index lexicographic prefix',
        'selected':history[:count],'unassessed':history[count:],'assessed_count':count,
        'certified_map_count':count,'assessed_unresolved_count':0,
        'unassessed_are_unresolved':False,'history_is_coverage_certificate':False}
    return closed,prior

class CoverageTests(unittest.TestCase):
    def test_exact_none_zero_one_and_unbounded_partitions(self):
        for limit in (None,0,1,100):
            a,b=fixture(limit);result=verify_history_coverage(a,b,limit,[])
            self.assertEqual(result['selected_count']+result['unassessed_count'],3)
    def test_shared_preclosure_state_and_metadata_are_exact(self):
        for location in ('round','state'):
            a,b=fixture(0)
            if location=='round':a['round']=4
            else:a['active_state']['next_needed']=[]
            with self.assertRaises(AssertionError):verify_history_coverage(a,b,0,[])
    def test_bad_partition_or_missing_coverage_rejected(self):
        for mutation in ('selection','missing','overlap'):
            a,b=fixture(1);c=a['active_state']['historical_candidate_collection']
            if mutation=='selection':c['selected']=[[1,0]]
            elif mutation=='overlap':c['unassessed'].insert(0,[0,1])
            else:del a['active_state']['historical_candidate_collection']
            with self.assertRaises((AssertionError,KeyError)):verify_history_coverage(a,b,1,[])
    def test_uncertified_status_and_invented_counts_rejected(self):
        for key,value in [('unassessed_are_unresolved',True),('history_is_coverage_certificate',True),('assessed_count',1),('certified_map_count',1)]:
            a,b=fixture(0);a['active_state']['historical_candidate_collection'][key]=value
            with self.assertRaises(AssertionError):verify_history_coverage(a,b,0,[])
    def test_retired_rows_must_be_selected_and_counted(self):
        a,b=fixture(1);c=a['active_state']['historical_candidate_collection'];c.update(certified_map_count=0,assessed_unresolved_count=1)
        row={'requested':[0,1],'unresolved':[{'integral':[2,0],'reason':'NoApplicableRule'}]}
        verify_history_coverage(a,b,1,[row])
        row['requested']=[1,0]
        with self.assertRaises(AssertionError):verify_history_coverage(a,b,1,[row])
    def test_option_and_sorted_unique_history_required(self):
        a,b=fixture(1)
        with self.assertRaises(AssertionError):verify_history_coverage(a,b,0,[])
        for row in (a,b):row['active_state']['historical_requests']=[[1,0],[0,1],[2,0]]
        with self.assertRaises(AssertionError):verify_history_coverage(a,b,1,[])

if __name__=='__main__':unittest.main()
