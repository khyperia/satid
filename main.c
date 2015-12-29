#include <predict/predict.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

const double DEG_TO_RAD = 0.01745329251;

typedef struct
{
    // observer
    double latitude;
    double longitude;
    double altitude;
    // target
    double rightAsc;
    double declination;
    // angle
    double maxAngle;
    // TLE file
    const char *tleFile;
    // time
    const char *time;
} arguments;

void usage()
{
    puts("Satid\n"
                 " -lat [degrees]  Latitude of observer\n"
                 " -lon [degrees]  Longitude of observer\n"
                 " -alt [meters]  Altitude of observer\n"
                 " -r [degrees]  Right ascension of satellite track\n"
                 " -d [degrees]  Declination of satellite track\n"
                 " -angle [degrees]  (optional, default: 1) Maximum angle at which results are no longer printed\n"
                 " -tle [filename]  Satellite info in the form of 3le (two-line elements with name as line 0)\n"
                 " -time [2000-01-01 12:00:00 -4]  End time of exposure (range is -10 +5 mins) (last num is timezone)\n"
                 "-r and -d can also be omitted, in this case, \"ra,dec\" number pairs are read from stdin\n"
                 "");
}

arguments parse_args(int argc, char **argv)
{
    arguments args;
    memset(&args, 0, sizeof(arguments));
    args.maxAngle = 1;
    for (int i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-lat"))
        {
            args.latitude = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "-lon"))
        {
            args.longitude = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "-alt"))
        {
            args.altitude = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "-r"))
        {
            args.rightAsc = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "-d"))
        {
            args.declination = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "-angle"))
        {
            args.maxAngle = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "-tle"))
        {
            args.tleFile = argv[++i];
        }
        else if (!strcmp(argv[i], "-time"))
        {
            args.time = argv[++i];
        }
        else
        {
            printf("Unknown argument %s, skipping\n", argv[i]);
        }
    }
    return args;
}

typedef struct
{
    predict_julian_date_t time_begin;
    predict_julian_date_t time_end;
    FILE *tleFile;
    predict_observer_t *observer;
    double maxAngle;
    double targetX, targetY, targetZ;
} satid_state;

satid_state *new_state(arguments args)
{
    satid_state result;
    struct tm parsedTime;
    int offTimezone;
    memset(&parsedTime, 0, sizeof(struct tm));
    if (sscanf(args.time, "%d-%d-%d %d:%d:%d %d",
               &parsedTime.tm_year, &parsedTime.tm_mon, &parsedTime.tm_mday,
               &parsedTime.tm_hour, &parsedTime.tm_min, &parsedTime.tm_sec,
               &offTimezone) != 7)
    {
        printf("Time was not in a proper format.\n");
        printf("Proper format is: yyyy-mm-dd hh:mm:ss tz\n");
        printf("(tz is an integer describing hour offset from local timezone)\n");
        return NULL;
    }
    parsedTime.tm_hour += offTimezone;
    parsedTime.tm_mon -= 1;
    parsedTime.tm_year -= 1900;
    predict_julian_date_t jul = predict_to_julian(mktime(&parsedTime));
    result.time_begin = jul + 5 / 1440.0;
    result.time_end = jul - 10 / 1440.0;
    result.tleFile = fopen(args.tleFile, "r");
    if (!result.tleFile)
    {
        printf("Could not open TLE file %s\n", args.tleFile);
        return NULL;
    }
    result.observer = predict_create_observer(
            "Observer", args.latitude * DEG_TO_RAD, args.longitude * DEG_TO_RAD, args.altitude);
    result.targetX = cos(args.declination * DEG_TO_RAD) * cos(args.rightAsc * DEG_TO_RAD);
    result.targetY = cos(args.declination * DEG_TO_RAD) * sin(args.rightAsc * DEG_TO_RAD);
    result.targetZ = sin(args.declination * DEG_TO_RAD);
    result.maxAngle = args.maxAngle;
    satid_state *ptr = malloc(sizeof(satid_state));
    if (!ptr)
    {
        printf("Out of memory\n");
        return NULL;
    }
    memcpy(ptr, &result, sizeof(satid_state));
    return ptr;
}

void delete_state(satid_state *state)
{
    if (state)
    {
        fclose(state->tleFile);
        predict_destroy_observer(state->observer);
        free(state);
    }
}

double search_at(satid_state *state, predict_orbit_t *orbit, predict_julian_date_t time)
{
    if (predict_orbit(orbit, time))
    {
        printf("predict_orbit returned nonzero\n");
        exit(1);
    }
    struct predict_observation observation;
    predict_observe_orbit(state->observer, orbit, &observation);
    observation.range_x /= observation.range;
    observation.range_y /= observation.range;
    observation.range_z /= observation.range;
    double cosAngle =
            observation.range_x * state->targetX +
            observation.range_y * state->targetY +
            observation.range_z * state->targetZ;
    return acos(cosAngle) / DEG_TO_RAD;
}

predict_julian_date_t optimize_time(satid_state *state, predict_orbit_t *orbit)
{
    predict_julian_date_t a = state->time_begin;
    predict_julian_date_t b = state->time_end;
    predict_julian_date_t bailCenter = (a + b) / 2;
    predict_julian_date_t bailRadius = a > b ? (a - b) / 2 : (b - a) / 2;
    a += (a - b) / 8; // expand range by a bit
    b += (b - a) / 8; // if it exits the range, break
    const predict_julian_date_t gr = (sqrt(5) - 1) / 2;
    predict_julian_date_t c = b - gr * (b - a);
    predict_julian_date_t d = a + gr * (b - a);
    double fc = search_at(state, orbit, c);
    double fd = search_at(state, orbit, d);
    while ((c > d ? c - d : d - c) > 0.00000001)
    {
        double currentBail = (a + b) / 2 - bailCenter;
        if ((currentBail < 0 ? -currentBail : currentBail) > bailRadius)
        {
            return 0;
        }
        if (fc < fd)
        {
            b = d;
            d = c;
            c = b - gr * (b - a);
            fd = fc;
            fc = search_at(state, orbit, c);
        }
        else
        {
            a = c;
            c = d;
            d = a + gr * (b - a);
            fc = fd;
            fd = search_at(state, orbit, d);
        }
    }
    return (b + a) / 2;
}

void strip_newline(char *str)
{
    size_t end = strlen(str) - 1;
    while (end > 0 && (str[end - 1] == '\r' || str[end - 1] == '\n'))
    {
        end--;
        str[end] = '\0';
    }
}

int perform_search(satid_state *state)
{
    char name[100];
    if (!fgets(name, sizeof(name), state->tleFile))
    {
        return 1;
    }
    strip_newline(name);
    char one[100];
    if (!fgets(one, sizeof(one), state->tleFile))
    {
        return 1;
    }
    strip_newline(name);
    char two[100];
    if (!fgets(two, sizeof(two), state->tleFile))
    {
        return 1;
    }
    strip_newline(name);
    char *tle[2] = {one, two};
    predict_orbit_t *orbit = predict_create_orbit(predict_parse_tle(tle));
    predict_julian_date_t time = optimize_time(state, orbit);
    if (time == 0)
    {
        return 0;
    }
    double angle = search_at(state, orbit, time);
    if (angle < state->maxAngle)
    {
        char timestr[50];
        time_t wholeTime = predict_from_julian(time);
        double fracTime = time - predict_to_julian(wholeTime);
        fracTime *= 86400;
        strftime(timestr, sizeof(timestr), "%Y-%m-%d %H:%M:%S", gmtime(&wholeTime));
        printf("Sat: %s (%d)\n", name, orbit->orbital_elements.satellite_number);
        printf("Angle: %f\n", angle);
        printf("Time: %s:%f\n", timestr, fracTime);
    }
    predict_destroy_orbit(orbit);
    return 0;
}

int main(int argc, char **argv)
{
    arguments args = parse_args(argc, argv);
    int notEnoughInfo = args.tleFile == NULL || args.time == NULL || args.latitude == 0 || args.longitude == 0 ||
                        args.altitude == 0;
    int noRaDec = args.rightAsc == 0 || args.declination == 0;
    if (notEnoughInfo)
    {
        usage();
        return 1;
    }
    if (noRaDec)
    {
        int count = 0;
        while (scanf("%lf,%lf", &args.rightAsc, &args.declination) == 2)
        {
            count++;
            printf(" -- RA/Dec -- %f,%f\n", args.rightAsc, args.declination);
            satid_state *state = new_state(args);
            if (!state)
            {
                return 1;
            }
            while (!perform_search(state))
            {
            }
            delete_state(state);
        }
        if (count == 0)
        {
            usage();
            return 1;
        }
        return 0;
    }
    else
    {
        satid_state *state = new_state(args);
        if (!state)
        {
            return 1;
        }
        while (!perform_search(state))
        {
        }
        delete_state(state);
        return 0;
    }
}