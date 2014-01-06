/*
 * Code section timer
 *
 *
 * */

#if !defined(TIMER_H)

  #define TIMER_H

  #if defined(TIMER)

    double timer_resolution(void);
    double timer_time(void);

    #if defined(PARALLEL)

      void timer_start(char *label);
      void timer_stop(char *label);

      #define TIMER_START(label) timer_start(label)
      #define TIMER_STOP(label) timer_stop(label)

    #else

      #define TIMER_START(label) fprintf(stderr, "timer start | %.9f | "label"\n", timer_time())
      #define TIMER_STOP(label) fprintf(stderr, "timer stop  | %.9f | "label"\n", timer_time())

    #endif

  #else

    #define TIMER_RESOLUTION
    #define TIMER_START
    #define TIMER_STOP

  #endif

#endif
