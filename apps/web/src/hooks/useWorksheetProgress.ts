/**
 * useWorksheetProgress - Hook for tracking worksheet completion progress.
 *
 * Persists progress state to localStorage so students can return to
 * where they left off. Tracks which sections are completed.
 *
 * @module hooks/useWorksheetProgress
 */

import { useState, useEffect, useCallback, useMemo } from 'react';

/**
 * Schema version for progress data.
 * Increment when making breaking changes to the schema.
 */
const PROGRESS_VERSION = 1;

/**
 * localStorage key prefix for progress data.
 */
const STORAGE_KEY_PREFIX = 'iqcp-worksheet-progress';

/**
 * Progress data stored in localStorage.
 */
export interface WorksheetProgress {
  /** Schema version for migrations */
  version: number;
  /** Lab pack identifier */
  labpack: string;
  /** Completion state for each section */
  sections: Record<
    string,
    {
      completed: boolean;
      completedAt?: string;
    }
  >;
  /** Last visited section ID */
  lastVisited?: string;
  /** Last update timestamp */
  lastUpdated: string;
}

/**
 * Return type for useWorksheetProgress hook.
 */
export interface UseWorksheetProgressReturn {
  /** Set of completed section IDs */
  completedSections: Set<string>;
  /** Check if a specific section is completed */
  isCompleted: (sectionId: string) => boolean;
  /** Mark a section as completed */
  markCompleted: (sectionId: string) => void;
  /** Mark a section as incomplete */
  markIncomplete: (sectionId: string) => void;
  /** Toggle section completion */
  toggleCompleted: (sectionId: string) => void;
  /** Clear all progress */
  clearProgress: () => void;
  /** Completion percentage (0-100) */
  progressPercent: number;
  /** Last visited section ID */
  lastVisited?: string;
  /** Update last visited section */
  setLastVisited: (sectionId: string) => void;
}

/**
 * Get storage key for a specific lab pack.
 */
function getStorageKey(labpack: string): string {
  return `${STORAGE_KEY_PREFIX}-${labpack}`;
}

/**
 * Load progress from localStorage.
 */
function loadProgress(labpack: string): WorksheetProgress | null {
  try {
    const key = getStorageKey(labpack);
    const stored = localStorage.getItem(key);

    if (!stored) return null;

    const data = JSON.parse(stored) as WorksheetProgress;

    // Check version and migrate if needed
    if (data.version !== PROGRESS_VERSION) {
      // For now, just clear old data on version mismatch
      localStorage.removeItem(key);
      return null;
    }

    return data;
  } catch (err) {
    console.error('Failed to load worksheet progress:', err);
    return null;
  }
}

/**
 * Save progress to localStorage.
 */
function saveProgress(labpack: string, progress: WorksheetProgress): void {
  try {
    const key = getStorageKey(labpack);
    localStorage.setItem(key, JSON.stringify(progress));
  } catch (err) {
    console.error('Failed to save worksheet progress:', err);
  }
}

/**
 * Create initial progress state.
 */
function createInitialProgress(labpack: string): WorksheetProgress {
  return {
    version: PROGRESS_VERSION,
    labpack,
    sections: {},
    lastUpdated: new Date().toISOString(),
  };
}

/**
 * useWorksheetProgress - Track worksheet completion progress.
 *
 * Persists to localStorage so progress survives page refreshes
 * and browser sessions. Provides methods to mark sections
 * complete/incomplete and calculate overall progress.
 *
 * @param labpack - Lab pack identifier (e.g., 'labpack1')
 * @param totalSections - Total number of trackable sections
 * @returns Progress state and methods
 *
 * @example
 * ```tsx
 * const { completedSections, toggleCompleted, progressPercent } =
 *   useWorksheetProgress('labpack1', 5);
 *
 * return (
 *   <div>
 *     <p>Progress: {progressPercent}%</p>
 *     <button onClick={() => toggleCompleted('section-1')}>
 *       {completedSections.has('section-1') ? 'Undo' : 'Complete'}
 *     </button>
 *   </div>
 * );
 * ```
 */
export function useWorksheetProgress(
  labpack: string,
  totalSections: number
): UseWorksheetProgressReturn {
  const [progress, setProgress] = useState<WorksheetProgress>(() => {
    return loadProgress(labpack) ?? createInitialProgress(labpack);
  });

  // Persist changes to localStorage
  useEffect(() => {
    saveProgress(labpack, progress);
  }, [labpack, progress]);

  // Compute completed sections as a Set
  const completedSections = useMemo(() => {
    return new Set(
      Object.entries(progress.sections)
        .filter(([, data]) => data.completed)
        .map(([id]) => id)
    );
  }, [progress.sections]);

  // Check if a section is completed
  const isCompleted = useCallback(
    (sectionId: string): boolean => {
      return progress.sections[sectionId]?.completed ?? false;
    },
    [progress.sections]
  );

  // Mark a section as completed
  const markCompleted = useCallback((sectionId: string) => {
    setProgress((prev) => ({
      ...prev,
      sections: {
        ...prev.sections,
        [sectionId]: {
          completed: true,
          completedAt: new Date().toISOString(),
        },
      },
      lastUpdated: new Date().toISOString(),
    }));
  }, []);

  // Mark a section as incomplete
  const markIncomplete = useCallback((sectionId: string) => {
    setProgress((prev) => ({
      ...prev,
      sections: {
        ...prev.sections,
        [sectionId]: {
          completed: false,
        },
      },
      lastUpdated: new Date().toISOString(),
    }));
  }, []);

  // Toggle section completion
  const toggleCompleted = useCallback(
    (sectionId: string) => {
      if (isCompleted(sectionId)) {
        markIncomplete(sectionId);
      } else {
        markCompleted(sectionId);
      }
    },
    [isCompleted, markCompleted, markIncomplete]
  );

  // Clear all progress
  const clearProgress = useCallback(() => {
    setProgress(createInitialProgress(labpack));
  }, [labpack]);

  // Calculate progress percentage
  const progressPercent = useMemo(() => {
    if (totalSections === 0) return 0;
    return Math.round((completedSections.size / totalSections) * 100);
  }, [completedSections.size, totalSections]);

  // Update last visited section
  const setLastVisited = useCallback((sectionId: string) => {
    setProgress((prev) => ({
      ...prev,
      lastVisited: sectionId,
      lastUpdated: new Date().toISOString(),
    }));
  }, []);

  return {
    completedSections,
    isCompleted,
    markCompleted,
    markIncomplete,
    toggleCompleted,
    clearProgress,
    progressPercent,
    lastVisited: progress.lastVisited,
    setLastVisited,
  };
}

export default useWorksheetProgress;
