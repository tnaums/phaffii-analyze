const std = @import("std");
const fasta = @import("root.zig");

pub fn main(init: std.process.Init) !void {
    const stdout = std.Io.File.stdout();

    const args = try init.minimal.args.toSlice(init.arena.allocator());
    if (args.len != 3) {
        std.debug.print("Usage: {s} <dna|protein> <filename>\n", .{args[0]});
        return;
    }

    var bmtype: fasta.biomolecule = undefined;
    if (std.meta.stringToEnum(fasta.biomolecule, args[1])) |value| {
        bmtype = value;
    } else {
        std.debug.print("No biomolecule type '{s}'\n", .{args[1]});
        return;
    }

    const filepath = args[2];
    const file = try std.Io.Dir.cwd().openFile(init.io, filepath, .{});
    defer file.close(init.io);

    var counter: u16 = 0;

    switch (bmtype) {
        .protein => {
            var queue: std.Io.Queue(fasta.Protein) = .init(&.{});
            var producer_task = try init.io.concurrent(fasta.parseProtein, .{ init.io, init.gpa, &queue, file });
            defer producer_task.cancel(init.io) catch {};
            var massiest: f32 = 0;

            const t_start = std.Io.Timestamp.now(init.io, .awake);
            while (true) {
                var myProtein = queue.getOne(init.io) catch |err| switch (err) {
                    error.Closed => break,
                    error.Canceled => return,
                };
                defer myProtein.deinit(init.gpa);
                counter += 1;

                if (myProtein.mass > massiest) {
                    massiest = myProtein.mass;
                    const printProtein = try std.fmt.allocPrint(init.gpa, "{f}\n", .{myProtein});
                    defer init.gpa.free(printProtein);
                    try stdout.writeStreamingAll(init.io, printProtein);
                }
            }
            const elapsed = t_start.durationTo(std.Io.Timestamp.now(init.io, .awake)).toMilliseconds();
            const elapsedPrint = try std.fmt.allocPrint(init.gpa, "elapsed time: {d} mS\n", .{elapsed});
            defer init.gpa.free(elapsedPrint);
            try stdout.writeStreamingAll(init.io, elapsedPrint);

            const massprint = try std.fmt.allocPrint(init.gpa, "Largest mass was: {d} kDa\n", .{massiest});
            defer init.gpa.free(massprint);
            try stdout.writeStreamingAll(init.io, massprint);
        },
        .dna => {
            // Starting DNA queue
            var queue: std.Io.Queue(fasta.DNA) = .init(&.{});
            var producer_task = try init.io.concurrent(fasta.parseDNA, .{ init.io, init.gpa, &queue, file });
            defer producer_task.cancel(init.io) catch {};

            var longest: u32 = 0;
            const t_start = std.Io.Timestamp.now(init.io, .awake);
            while (true) {
                var myDNA = queue.getOne(init.io) catch |err| switch (err) {
                    error.Closed => break,
                    error.Canceled => return,
                };
                defer myDNA.deinit(init.gpa);
                counter += 1;

                const length: u32 = @intCast(myDNA.sequence.len);
                if (length > longest) {
                    longest = length;
                }

                const fs = try std.fmt.allocPrint(init.gpa, "{f}\n", .{myDNA});
                defer init.gpa.free(fs);
                try stdout.writeStreamingAll(init.io, fs);

                try myDNA.addTranslation(init.gpa);
                myDNA.printOrfs(init.gpa, 50);
                try myDNA.mapDNA(init.gpa, init.io, stdout, fasta.readingframe.second);

                try stdout.writeStreamingAll(init.io, "------------------------------------------------------------\n");
            }
            const elapsed = t_start.durationTo(std.Io.Timestamp.now(init.io, .awake)).toMilliseconds();
            const elapsedPrint = try std.fmt.allocPrint(init.gpa, "elapsed time: {d} mS\n", .{elapsed});
            defer init.gpa.free(elapsedPrint);
            try stdout.writeStreamingAll(init.io, elapsedPrint);

            const longprint = try std.fmt.allocPrint(init.gpa, "Longest gene was {d} nucleotides\n", .{longest});
            defer init.gpa.free(longprint);
            try stdout.writeStreamingAll(init.io, longprint);
        },
    }

    const finalTally = try std.fmt.allocPrint(init.gpa, "Created {d} Fasta objects\n", .{counter});
    defer init.gpa.free(finalTally);
    try stdout.writeStreamingAll(init.io, finalTally);


}
