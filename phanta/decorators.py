import warnings
import functools


def experimental(description=None):
    """
    MARKS FUNCTION AS EXPERIMENTAL:
        - Raises a warning to indicate that this function is experimental.
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            message = f"The function '{func.__name__}' is experimental and has not been tested."
            if description:
                message += f" {description}"
            warnings.warn(
                message,
                category=UserWarning,
                stacklevel=2
            )
            return func(*args, **kwargs)
        return wrapper
    return decorator


def non_gmp_compliant(description=None):
    """
    MARKS FUNCTION:
        - Raises a warning to indicate that this function does not comply with Good Automated Manufacturing Practices
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            message = f"The function '{func.__name__}' is not GAMP compliant."
            if description:
                message += f" {description}"
            warnings.warn(
                message,
                category=UserWarning,
                stacklevel=2
            )
            return func(*args, **kwargs)
        return wrapper
    return decorator
