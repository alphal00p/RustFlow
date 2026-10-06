"""Run actual Future polling against old and fixed helper ASTs without workers."""
import ast
from concurrent.futures import Future, TimeoutError as FutureTimeoutError
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

path = Path('examples/hep/gg_hg_acceptance.py')
old = subprocess.check_output(['git', 'show', '3759d72f:' + str(path)], text=True)
current = path.read_text()

def load_wait(source):
    node = next(node for node in ast.parse(source).body if isinstance(node, ast.FunctionDef) and node.name == 'wait_for_stage')
    namespace = {'FutureTimeoutError': FutureTimeoutError}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), 'exec'), namespace)
    return namespace['wait_for_stage']

def pending(wait_for_stage):
    future = Future()
    value = object()
    updates = []
    def publish(state):
        updates.append(state['done'])
        if len(updates) == 2:
            future.set_result(value)
    session = SimpleNamespace(snapshot=lambda: {'done': future.done()}, wait=future.result)
    try:
        result = wait_for_stage(session, publish, poll_interval=0.001)
    except FutureTimeoutError:
        return {'result': 'poll timeout escaped', 'progress': updates}
    assert result is value
    return {'result': 'completed', 'progress': updates}

def worker_timeout(wait_for_stage):
    future = Future()
    original = TimeoutError('worker failure')
    updates = []
    def wait(timeout=None):
        if not future.done():
            future.set_exception(original)
        return future.result(timeout)
    session = SimpleNamespace(snapshot=lambda: {'done': future.done()}, wait=wait)
    try:
        wait_for_stage(session, lambda state: updates.append(state['done']), poll_interval=0.001)
    except TimeoutError as error:
        assert error is original
    else:
        raise AssertionError('worker timeout swallowed')
    return updates

results = {
    'python': sys.version,
    'timeout_classes_are_aliases': FutureTimeoutError is TimeoutError,
    'old_pending': pending(load_wait(old)),
    'fixed_pending': pending(load_wait(current)),
    'old_worker_failure_progress': worker_timeout(load_wait(old)),
    'fixed_worker_failure_progress': worker_timeout(load_wait(current)),
}
assert results['old_pending']['result'] == ('completed' if FutureTimeoutError is TimeoutError else 'poll timeout escaped')
assert results['fixed_pending'] == {'result': 'completed', 'progress': [False, False, True]}
assert results['old_worker_failure_progress'] == [False]
assert results['fixed_worker_failure_progress'] == [False, True]
print(json.dumps(results, indent=2))
