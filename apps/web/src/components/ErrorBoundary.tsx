/**
 * ErrorBoundary - React error boundary component.
 *
 * Catches JavaScript errors in child component tree and displays
 * a fallback UI instead of crashing the whole app.
 *
 * @module components/ErrorBoundary
 */

import { Component, type ReactNode, type ErrorInfo } from 'react';

interface ErrorBoundaryProps {
  children: ReactNode;
  /** Optional fallback component */
  fallback?: ReactNode;
}

interface ErrorBoundaryState {
  hasError: boolean;
  error: Error | null;
  errorInfo: ErrorInfo | null;
}

/**
 * Error boundary component for catching React errors.
 *
 * @example
 * ```tsx
 * <ErrorBoundary>
 *   <MyComponent />
 * </ErrorBoundary>
 * ```
 */
export class ErrorBoundary extends Component<ErrorBoundaryProps, ErrorBoundaryState> {
  constructor(props: ErrorBoundaryProps) {
    super(props);
    this.state = { hasError: false, error: null, errorInfo: null };
  }

  static getDerivedStateFromError(error: Error): Partial<ErrorBoundaryState> {
    return { hasError: true, error };
  }

  componentDidCatch(error: Error, errorInfo: ErrorInfo): void {
    console.error('ErrorBoundary caught an error:', error, errorInfo);
    this.setState({ errorInfo });
  }

  render(): ReactNode {
    if (this.state.hasError) {
      if (this.props.fallback) {
        return this.props.fallback;
      }

      return (
        <div className="bg-red-50 border border-red-200 rounded-xl p-6 m-4">
          <h2 className="text-lg font-semibold text-red-800 mb-2">
            Something went wrong
          </h2>
          <p className="text-red-600 mb-4">
            An error occurred while rendering this component.
          </p>
          {this.state.error && (
            <div className="bg-red-100 rounded-lg p-4 mb-4">
              <p className="font-mono text-sm text-red-700 break-all">
                {this.state.error.toString()}
              </p>
            </div>
          )}
          {this.state.errorInfo && (
            <details className="mt-4">
              <summary className="cursor-pointer text-sm text-red-600 hover:text-red-800">
                Component Stack
              </summary>
              <pre className="mt-2 text-xs text-red-600 overflow-auto max-h-48 bg-red-100 p-2 rounded">
                {this.state.errorInfo.componentStack}
              </pre>
            </details>
          )}
          <button
            onClick={() => window.location.reload()}
            className="mt-4 px-4 py-2 bg-red-600 text-white rounded-lg hover:bg-red-700 focus:outline-none focus:ring-2 focus:ring-red-500"
          >
            Reload Page
          </button>
        </div>
      );
    }

    return this.props.children;
  }
}

export default ErrorBoundary;
