def test(func):
    """Marks a function as a test case.
    
    Example:

    @minitest.test
    def test1():
        assert True
    """
    test_list = func.__globals__.setdefault('_minitests', [])
    test_list.append(func)
    return func

def interactive(func):
    """Marks a function as an interactive test case.
    
    Example:

    @minitest.interactive
    def test2():
        answer = input('Type yes: ')
        assert answer == 'yes'
    """
    test_list = func.__globals__.setdefault('_minitests_i', [])
    test_list.append(func)
    return func

def runner(func):
    """Turns a function into a test runner.

    The function must take a single argument (controlling whether interactive
    tests are executed). The body of the function will be executed first,
    followed by all the declared test cases in the current module.
    """
    test_list = func.__globals__.setdefault('_minitests', [])
    i_test_list = func.__globals__.setdefault('_minitests_i', [])
    module_name = func.__module__
    def new_func(interactive):
        ret = func(interactive)
        print('Running tests in', module_name)
        for f in test_list:
            f()
        if interactive:
            for f in i_test_list:
                f()
        return ret
    return new_func
