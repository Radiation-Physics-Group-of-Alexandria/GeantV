//===--- GeantTrackStat.h - Geant-V -----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantTrackStat.h
 * @brief Implementation of statistics for track array in Geant-V prototype
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_TRACK_STAT
#define GEANT_TRACK_STAT

#include <mutex>

#include "GeantFwd.h"

/**
 * @brief Class GeantTrackStat
 * @details Used statistic object for a track array
 *
 */
class GeantTrackStat {
public:
  using GeantTrack = Geant::GeantTrack;
  using GeantTrack_v = Geant::GeantTrack_v;
  using GeantTaskData = Geant::GeantTaskData;

  int fNslots;   /** Number of event slots */
  int *fNtracks; /** [fNslots] Number of tracks from an event */
  int *fNsteps;  /**[fNslots] Cumulated number of steps per event */
  std::mutex fMutex; /** Mutex */

private:
  /**
   * @brief Copy constructor
   */
  GeantTrackStat(const GeantTrackStat &other);

  /**
   * @brief Operator =
   *
   * @param other ?????
   */
  GeantTrackStat &operator=(const GeantTrackStat &other);

public:
  /** @brief GeantTrackStat constructor */
  GeantTrackStat() = default;

  /**
   * @brief GeantTrackStat constructor
   *
   * @param nslots Number of event slots
   */
  GeantTrackStat(int nslots);

  /** @brief GeantTrackStat destructor */
  virtual ~GeantTrackStat();

  /**
   * @brief Operator +=
   *
   * @param other ??????
   */
  GeantTrackStat &operator+=(const GeantTrackStat &other);

  /**
   * @brief Operator -=
   *
   * @param other ?????
   */
  GeantTrackStat &operator-=(const GeantTrackStat &other);

  /**
   * @brief Function that add track
   *
   * @param track Track that should be added
   */
  void AddTrack(const GeantTrack &track);

  /**
   * @brief Function that add track from track_v array
   *
   * @param trackv Track that should be added
   * @param itr Track ID ?????
   */
  void AddTrack(const GeantTrack_v &trackv, int itr);

  /**
   * @brief Function that add tracks from track_v array
   *
   * @param trackv Tracks that should be added
   */
  void AddTracks(const GeantTrack_v &trackv);

  /**
   * @brief Function that remove tracks
   *
   * @param trackv Tracks that should be deleted
   */
  void RemoveTracks(const GeantTrack_v &trackv);

  /**
   * @brief Function that init array of event slots ?????
   *
   * @param nslots Number of event slots
   */
  void InitArrays(int nslots);

  /** @brief Print function */
  void Print(const char *option = "") const;

  /** @brief Reset function */
  void Reset();

};
#endif
