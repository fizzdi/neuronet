#pragma once

#define MAX_STR_LEN 5000
// Максимальная длина сообщения лога
#define LOG_MSG_LEN 100000
// Идентификаторы исключений
/*******************************************/
#define E_COMPLETION_PORT_ASSOCIATION  1
#define E_READ_JOB_MEMORY_INFORMATION  2
#define E_TERMINATE_JOB_OBJECT         3
#define E_CREATE_JOB_OBJECT            4
#define E_CREATE_COMPLETION_PORT       5
#define E_SET_LIMITS                   6
#define E_ASSIGN_PROCESS_TO_JOB        7
#define E_CLOSE_ALL_FOR_JOB            8
#define E_MEMORY_LIMIT                 9
#define E_RESUME_THREAD               10
#define E_READ_COMPLETION_PORT        11
#define E_CREATE_PROCESS              12
#define E_SUCCESSFULL                 13
#define E_COMPILE_ERROR               14
#define E_SET_JOB_INFORMATION         15
#define E_NO_TESTS                    16
#define E_COPY_INPUT_FILE             17
#define E_LIBRARY_ACCESS              18
#define E_READ_DLL_FUNCTION           19
#define E_CREATE_PIPE                 20
#define E_CREATE_COMPILATION_PROCESS  21
#define E_SET_JOB_LIMITS              22
#define E_SET_JOB_RESTRICTIONS        23
#define E_SET_JOB_TIME_INFORMATION    24
#define E_LOCK_TABLE_SQL              25
#define E_SELECT_SOLUTION_SQL         26
#define E_NO_SOLUTION                 27
#define E_UPDATE_CHECK_STATE          28
#define E_CREATE_THREAD               29
#define E_INIT_TESTER                 30
#define E_PASSWORD_REQUEST            31
#define E_INTERSEPT_API               32
#define E_TIME_LIMIT                  33
#define E_RUNTIME_ERROR				  34
/*******************************************/

// Количество исключений
#define EXCEPTION_COUNT 40

// Константа - массив сообщений для описанных выше исключений
/*****************************************************************************************************/
const char Messages[EXCEPTION_COUNT][MAX_STR_LEN] = { "Error associating completion port with job",
"Error reading job memory information",
"Error terminating job object",
"Error creating job object",
"Error creating completion port",
"",
"Error assigning process to job object",
"",
"",
"Error resuming thread",
"Error reading completion port informaition",
"Creating process error",
"",
"",
"",
"",
"Error coping input file",
"Library access error",
"Error reading dll function",
"Creating pipe error",
"Creating compilation process error",
"Error setting limit information for job",
"Error setting job basic restrictions",
"Error setting end of job time information",
"Lock table SQL query error",
"Select solution SQL query error",
"",
"Update check state SQL query error",
"Can't create working thread",
"",
"",
"Interception API error",
"" };
/*****************************************************************************************************/
