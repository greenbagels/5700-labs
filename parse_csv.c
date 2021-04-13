
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct bitscope_data
{
    int channels;
    size_t samples;
    size_t count;
    double rate;
    double *data;
};

int parse_csv(FILE*, struct bitscope_data*);
double *gen_time_array(struct bitscope_data*);
int calc_avgs(struct bitscope_data*, double*[2]);
int print_chans(char*, int, size_t, double*, double*[2]);

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        fprintf(stderr, "Error: expected a filename argument");
        return 1;
    }
    
    // Just for generality sake
    char filename[80];
    snprintf(filename, 80, "%s%s", argv[1], ".csv"); 
    FILE* rlc_input = fopen(filename, "r");
    if (rlc_input == NULL)
    {
        exit(3);
        // TODO: error handling
    }
    
    struct bitscope_data bsdata;
    // Parse the oscilloscope data into the struct
    parse_csv(rlc_input, &bsdata);
    // We're done parsing the CSV, close up shop
    fclose(rlc_input);

    // Generate the list of times using the rate and count
    double* times = gen_time_array(&bsdata);
    if (times == NULL)
    {
        // TODO: error handling
    }
    // TODO: abstract this
    double* chans[2];
    for (int i = 0; i < 2; i++)
    {
        chans[i] = (i < bsdata.channels ? (double*)malloc(bsdata.count * 2 * sizeof(double)) : NULL);
    }

    calc_avgs(&bsdata, chans);
    print_chans(argv[1], bsdata.channels, bsdata.count, times, chans);
}

/*! Gets the number of lines in the current file.
 *
 *  This function is guaranteed to leave the file position indicator unmodified.
 */
int get_file_len(FILE* stream, size_t *flen)
{
    // First, record the current position of the stream.
    long pos = ftell(stream);
    // Now, rewind to the beginning and parse.
    rewind(stream);

    char* temp = NULL;
    size_t line_len = 0u;
    size_t len = 0u;

    while (getline(&temp, &line_len, stream) > 0)
    {
        len++;
    }

    // Cleanup what getline allocated.
    free(temp);
    // Leave the stream how it was before we started working on it
    fseek(stream, pos, SEEK_SET);

    *flen = len;

    return 0;
}

int get_num_channels(FILE* stream, int* num_channels)
{
    // First, record the current position of the stream.
    long pos = ftell(stream);
    // Now, rewind to the beginning and parse.
    rewind(stream);

    size_t len = 0;
    char *str = NULL;
    // Skip first line, which just has labels
    getline(&str, &len, stream);

    // Check the Channel indices of the 2nd and third lines.
    int channels[2];
    for (int i = 0; i < 2; i++)
    {
        getline(&str, &len, stream);
        // Skip the TRIGGER and STAMP fields
        strtok(str, ",");
        strtok(NULL, ",");
        // Now parse the CHANNEL field
        channels[i] = atoi(strtok(NULL, ","));
    }

    free(str);
    fseek(stream, pos, SEEK_SET);

    *num_channels = (channels[1] == channels[0] ? 1 : 2);
    return 0;
}

int get_rate_count(FILE* stream, double *rate, size_t *count)
{
    long pos = ftell(stream);
    rewind(stream);

    size_t len = 0;
    char *str = NULL;

    getline(&str, &len, stream);
    getline(&str, &len, stream);

    strtok(str, ",");  // TRIGGER
    strtok(NULL, ","); // STAMP
    strtok(NULL, ","); // CHANNEL
    strtok(NULL, ","); // INDEX
    strtok(NULL, ","); // TYPE
    strtok(NULL, ","); // DELAY
    strtok(NULL, ","); // FACTOR

    *rate = strtod(strtok(NULL, ","), NULL);
    *count = strtol(strtok(NULL, ","), NULL, 10);

    free(str);
    fseek(stream, pos, SEEK_SET);

    return 0;
}

int parse_csv(FILE* stream, struct bitscope_data* data)
{
    // Some useful metadata for later...
    get_num_channels(stream, &data->channels);
    size_t lines;
    get_file_len(stream, &lines);
    data->samples = (lines-1)/data->channels;
    get_rate_count(stream, &data->rate, &data->count);
    // Supposedly the last few data points are bad, so we should ignore them
    data->count -= 5;

    fprintf(stderr, "%d %lu %lu %lf\n", data->channels, data->samples, data->count, data->rate);

    data->data = (double*)malloc(data->samples * data->count * data->channels * sizeof(double));
    fprintf(stderr, "data address :%p\n", data->data);

    // Now, loop to read in the data.

    size_t len = 0;
    char *str = NULL;

    // Skip the first line, which only has labels
    getline(&str, &len, stream);

    for (size_t i = 0; i < data->samples; i++)
    {
        // Get every line, which alternates channels before samples
        getline(&str, &len, stream);

        // Discard the first 2 fields (TRIGGER thru STAMP)
        strtok(str, ",");
        strtok(NULL, ",");

        // Parse the next field (CHANNEL)
        size_t j = strtol(strtok(NULL, ","), NULL, 10);

        // Discard the next 6 fields (INDEX thru COUNT)
        for (int k = 0; k < 6; k++)
        {
            strtok(NULL, ",");
        }

        // printf("%u %u\n", i, j);

        // Parse the last field (DATA)
        for (size_t k = 0; k < data->count; k++)
        {
            data->data[j*data->count*data->samples + k*data->samples + i] = strtod(strtok(NULL, ",\r\n"), NULL);
        }
    }
    return 0;
}

double* gen_time_array(struct bitscope_data* data)
{
    double *times = (double*)malloc((data->count) * sizeof(double));
    for (size_t i = 0; i < data->count; i++)
    {
        times[i] = i / data->rate;
    }

    return times;
}

int calc_avgs(struct bitscope_data* data, double *chans[2])
{

    fprintf(stderr, "%d %lu %lu %lf\n", data->channels, data->samples, data->count, data->rate);
    fprintf(stderr, "data address :%p\n", data->data);

    for (int i = 0; i < data->channels; i++)
    {
        for (size_t j = 0; j < data->count; j++)
        {
            double avg = 0.;
            double err = 0.;
            for (size_t k = 0; k < data->samples; k++)
            {
                avg += data->data[i*data->samples * data->count + j*data->samples + k];
            }
            avg /= data->samples;
            for (size_t k = 0; k < data->samples; k++)
            {
                err += pow(avg - data->data[i*data->samples*data->count + j*data->samples + k], 2);
            }
            err = sqrt(err/data->samples/data->samples);
            chans[i][j] = avg;
            chans[i][data->count + j] = err;
        }
    }
    return 0;
}

int print_chans(char *prefix, int numchans, size_t count, double* times, double* chans[2])
{

    char fname[20];
    snprintf(fname, sizeof(fname), "%s.dat", prefix);
    FILE* ofile = fopen(fname, "w");
    if (ofile == NULL)
    {
        // TODO: better error handling
        exit(2);
    }
    for (size_t i = 0; i < count; i++)
    {
        for (size_t j = 0; j < numchans; j++)
        {
            // Fails on the loop i=j=0
            // Also yes i know the loop ordering is suboptimal
            if (fprintf(ofile, "%lf %lf %lf ", times[i], chans[j][i], chans[j][i+count]) < 0)
            {
                // TODO: better error handling
                exit(1);
            }
        }
        fprintf(ofile, "\n");
    }
    fclose(ofile);
    return 0;
}
