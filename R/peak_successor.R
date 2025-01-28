
# Finds which labeled peaks are in same location as peak in previous curve
peak_successor <- function(peaks1, valleys1, peaks2, valleys2, labels1, t_grid) {
  # Helper function to compute ranges for peaks
  compute_peak_ranges <- function(peaks, valleys) {
    # Ensure valleys include endpoints
    all_valleys <- sort(unique(c(t_grid[1], valleys, t_grid[length(t_grid)])))

    # For each peak, find its containing range
    ranges <- lapply(peaks, function(peak) {
      left_valley <- max(all_valleys[all_valleys < peak])
      right_valley <- min(all_valleys[all_valleys > peak])
      c(left_valley, right_valley)
    })

    # Convert list to matrix
    do.call(rbind, ranges)
  }

  # Helper function to assign labels to peaks based on ranges
  assign_labels_to_peaks <- function(peaks2, ranges, peaks1, labels1) {
    labels <- numeric(length(peaks2))

    for(i in seq_along(peaks2)) {
      # Find which ranges contain this peak
      in_range <- peaks2[i] >= ranges[,1] & peaks2[i] <= ranges[,2]
      matching_ranges <- which(in_range)

      if(length(matching_ranges) > 0) {
        if(length(matching_ranges) > 1) {
          # If multiple ranges match, choose the closest peak
          distances <- abs(peaks1[matching_ranges] - peaks2[i])
          closest_idx <- matching_ranges[which.min(distances)]
          labels[i] <- labels1[closest_idx]
        } else {
          labels[i] <- labels1[matching_ranges]
        }
      }
    }
    labels
  }

  # Helper function to resolve overlapping labels
  resolve_overlapping_labels <- function(labels, peaks1, peaks2, labels1) {
    unique_labels <- unique(labels[labels > 0])

    for(label in unique_labels) {
      duplicates <- which(labels == label)
      if(length(duplicates) > 1) {
        # Keep the closest peak and reset others
        orig_peak <- peaks1[which(labels1 == label)]
        distances <- abs(orig_peak - peaks2[duplicates])
        keep_idx <- duplicates[which.min(distances)]
        remove_idx <- duplicates[duplicates != keep_idx]
        labels[remove_idx] <- 0
      }
    }
    labels
  }

  # Main function logic
  if(length(peaks1) == 0) {
    # If no peaks in first function, assign new labels to all peaks in second
    labelMax <- max(labels1, 0)
    labels2 <- labelMax + seq_along(peaks2)
    return(list(
      labels = labels2,
      labelMax = labelMax + length(peaks2)
    ))
  }

  # Compute ranges for peaks in first function
  ranges <- compute_peak_ranges(peaks1, valleys1)

  # Assign initial labels to peaks in second function
  labels2 <- assign_labels_to_peaks(peaks2, ranges, peaks1, labels1)

  # Resolve any overlapping labels
  labels2 <- resolve_overlapping_labels(labels2, peaks1, peaks2, labels1)

  # Assign new labels to unmatched peaks
  labelMax <- max(labels1)
  unmatched <- which(labels2 == 0)
  if(length(unmatched) > 0) {
    labels2[unmatched] <- labelMax + seq_along(unmatched)
    labelMax <- labelMax + length(unmatched)
  }

  # Return both new labels and updated labelMax
  return(list(
    labels = labels2,
    labelMax = labelMax
  ))
}
