const std = @import("std");
const fasta = @import("fasta");
const Io = std.Io;

const promoter = enum { aox1, gap, unknown };
const secretion = enum { alpha, ost, cytoplasmic };
const organism = enum { unknown, verticillioides };

const KpError = error{
    PromoterNotFound,
    CodingNotFound,
};

const BestLocalBlastHit = struct {
    identity: usize,
    def: []const u8,

    pub fn init(allocator: std.mem.Allocator, idLength: usize, def: []const u8) !BestLocalBlastHit {
        const d = try allocator.dupe(u8, def);
        errdefer allocator.free(d);
        return .{
            .identity = idLength,
            .def = d,
        };
    }

    pub fn deinit(self: BestLocalBlastHit, allocator: std.mem.Allocator) void {
        allocator.free(self.def);
    }
};

const Komagataella = struct {
    plasmid: fasta.DNA,
    promoter: promoter,
    coding: fasta.DNA,
    secretion: secretion,
    recombinant: fasta.Protein,
    localBlast: ?BestLocalBlastHit,

    pub fn init(io: Io, allocator: std.mem.Allocator, filepath: []const u8) !Komagataella {
        var dna = try fasta.parseSingleDNA(io, allocator, filepath);
        errdefer dna.deinit(allocator);
        try dna.addTranslation(allocator);
        const p = checkPromoter(dna);
        var codingRegion = try getCodingRegion(allocator, dna, p);
        errdefer codingRegion.deinit(allocator);
        try codingRegion.addTranslation(allocator);
        const s = determineSecretion(codingRegion);
        const r = try recombinantProtein(allocator, codingRegion, s);
        errdefer r.deinit(allocator);
        const lb = try localBlastBestHit(io, allocator, r);

        return Komagataella{
            .plasmid = dna,
            .promoter = p,
            .coding = codingRegion,
            .secretion = s,
            .recombinant = r,
            .localBlast = lb,
        };
    }

    pub fn deinit(self: Komagataella, allocator: std.mem.Allocator) void {
        self.plasmid.deinit(allocator);
        self.coding.deinit(allocator);
        self.recombinant.deinit(allocator);
        if (self.localBlast) |value| {
            value.deinit(allocator);
        }
    }

    fn localBlastBestHit(io: Io, allocator: std.mem.Allocator, r: fasta.Protein) !?BestLocalBlastHit {
        // first, write mature recombinant protein to disk as a fasta file
        // this will be overwritten each time
        const queryfile = try std.Io.Dir.cwd().createFile(io, "blastdb/query.fa", .{});
        defer queryfile.close(io);

        var writebuf: [1024]u8 = undefined;
        const queryPrint = try std.fmt.bufPrint(&writebuf, "{f}", .{r});
        try queryfile.writeStreamingAll(io, queryPrint);

        // second, perform localblast against database
        const argv1 = [_][]const u8{ "blastp", "-query", "blastdb/query.fa", "-db", "blastdb/combined/combined.aa.fasta", "-outfmt", "5", "-out", "-", "-evalue", "1e-50" };
        const proc1 = try std.process.run(
            allocator,
            io,
            .{ .argv = &argv1 },
        );
        defer allocator.free(proc1.stdout);
        defer allocator.free(proc1.stderr);
        switch (proc1.term) {
            .exited => |eCode| {
                if (eCode != 0) {
                    std.debug.print(
                        "P1 did not exit cleanly with code {}\n",
                        .{eCode},
                    );
                } else {
                    const localBlast = try getBlastpBestHit(allocator, proc1.stdout);
                    if (localBlast) |value| {
                        return value;
                    }
                    return null;
                }
            },
            .signal => {
                std.debug.print("P1 was killed by a signal.\n", .{});
            },
            .stopped => {
                std.debug.print("P1 was stopped.\n", .{});
            },
            .unknown => {
                std.debug.print("P1 terminated: unknown status.\n", .{});
            },
        }
        return null;
    }

    fn checkPromoter(dna: fasta.DNA) promoter {
        const aox1Start = std.mem.indexOf(u8, dna.sequence, "AGATCTAACATCCAAA");
        const aox1End = std.mem.indexOf(u8, dna.sequence, "ATTCGAAACGA") orelse 0;

        const gapStart = std.mem.indexOf(u8, dna.sequence, "AGATCTTTTTTGTAGAAATG");
        const gapEnd = std.mem.indexOf(u8, dna.sequence, "TATTTCAATCAATTGAACAAC") orelse 0;

        if (aox1Start) |value| {
            if (value == 0 and aox1End == 930) {
                return promoter.aox1;
            }
        }
        if (gapStart) |value| {
            if (value == 0 and gapEnd == 459) {
                return promoter.gap;
            }
        }
        return promoter.unknown;
    }

    fn getCodingRegion(allocator: std.mem.Allocator, dna: fasta.DNA, p: promoter) !fasta.DNA {
        const newHeader = try std.fmt.allocPrint(allocator, "{s}|{s}", .{ dna.header, "coding" });
        defer allocator.free(newHeader);

        switch (p) {
            .aox1 => {
                const codingEnd = std.mem.indexOf(u8, dna.sequence, "GTTTGTAGCCTTAGA") orelse dna.sequence.len;
                return try fasta.DNA.init(allocator, newHeader, dna.sequence[940..codingEnd]);
            },
            .gap => {
                const codingEnd = std.mem.indexOf(u8, dna.sequence, "GTTTTAGCCTTAGAC") orelse dna.sequence.len;
                return try fasta.DNA.init(allocator, newHeader, dna.sequence[492..codingEnd]);
            },
            .unknown => {
                return KpError.PromoterNotFound;
            },
        }
    }

    fn determineSecretion(dna: fasta.DNA) secretion {
        if (std.mem.startsWith(u8, dna.sequence, "ATGAGATTTCCTTCA")) {
            return secretion.alpha;
        }
        if (std.mem.startsWith(u8, dna.sequence, "ATGAGGCAGGTTTGG")) {
            return secretion.ost;
        }
        return secretion.cytoplasmic;
    }

    fn recombinantProtein(allocator: std.mem.Allocator, codingRegion: fasta.DNA, s: secretion) !fasta.Protein {
        const newHeader = try std.fmt.allocPrint(allocator, "{s}|{s}", .{ codingRegion.header[0 .. codingRegion.header.len - 7], "mature" });
        defer allocator.free(newHeader);
        var fullP: []u8 = undefined;

        if (codingRegion.translation) |value| {
            fullP = value[0];
        }
        var parts = std.mem.tokenizeScalar(u8, fullP, '*');
        var shortLeft = parts.next() orelse fullP;

        switch (s) {
            .alpha => {
                return try fasta.Protein.init(allocator, newHeader, shortLeft[89..]);
            },

            .ost => {
                return try fasta.Protein.init(allocator, newHeader, shortLeft[92..]);
            },
            .cytoplasmic => {
                return try fasta.Protein.init(allocator, newHeader, shortLeft);
            },
        }
    }

    pub fn format(self: Komagataella, writer: *Io.Writer) !void {
        const bp = self.plasmid.sequence.len;
        try writer.print("   plasmid: {s}\n    length: {d}bp\n", .{ self.plasmid.header, bp });
        switch (self.promoter) {
            .aox1 => {
                try writer.print("  promoter: methanol inducible aox1\n", .{});
            },
            .gap => {
                try writer.print("  promoter: constitutive gap promoter\n", .{});
            },
            .unknown => {
                try writer.print("  promoter: unknown, probably not a pPICZ expression plasmid\n", .{});
            },
        }
        switch (self.secretion) {
            .alpha => {
                try writer.print(" secretion: α factor\n", .{});
            },
            .ost => {
                try writer.print(" secretion: ost1\n", .{});
            },
            .cytoplasmic => {
                try writer.print(" secretion: standard SSS not detected, may be cytoplasmic\n", .{});
            },
        }
        if (self.localBlast) |value| {
            try writer.print("localBlast: {s}\n  identity: {d}/{d}\n", .{ value.def, value.identity, self.recombinant.sequence.len });
        }
        try writer.print("\n", .{});
        try writer.print("{f}\n", .{self.recombinant});
    }

    pub fn printCSV(self: Komagataella, io: Io, file: Io.File) !void {
        var buffer: [70]u8 = undefined;
        const line = try std.fmt.bufPrint(&buffer, "{s},{s},{s},{d:.2},{d}\n", .{ self.plasmid.header[0..7], @tagName(self.promoter), @tagName(self.secretion), self.recombinant.mass, self.recombinant.sequence.len });
        try file.writeStreamingAll(io, line);
    }
};

pub fn main(init: std.process.Init) !void {
    const stdout = std.Io.File.stdout();

    const args = try init.minimal.args.toSlice(init.arena.allocator());
    if (args.len != 2) {
        std.debug.print("Usage: {s} <filename.fa>\n", .{args[0]});
        return;
    }

    const filepath = args[1];

    var k = try Komagataella.init(init.io, init.gpa, filepath);
    defer k.deinit(init.gpa);

    const fs = try std.fmt.allocPrint(init.gpa, "{f}\n", .{k.coding});
    defer init.gpa.free(fs);
    try stdout.writeStreamingAll(init.io, fs);
        
    try k.coding.mapDNA(init.gpa, init.io, stdout, fasta.readingframe.first);

    try stdout.writeStreamingAll(init.io, "----------------------------------------\n");

    const kf = try std.fmt.allocPrint(init.gpa, "{f}\n", .{k});
    defer init.gpa.free(kf);
    try stdout.writeStreamingAll(init.io, kf);

    try k.printCSV(init.io, stdout);
}

fn getBlastpBestHit(allocator: std.mem.Allocator, blastResult: []u8) !?BestLocalBlastHit {
    var def: []u8 = undefined;
    var size: usize = undefined;

    //    std.debug.print("{s}\n", .{blastResult});
    var hits = std.mem.tokenizeSequence(u8, blastResult, "<Hit>");
    _ = hits.next(); // discard header info
    if (hits.next()) |value| { // if there is a hit, process the best (first) one
        var lines = std.mem.tokenizeScalar(u8, value, '\n');
        while (lines.next()) |line| {
            const cleaned = std.mem.trim(u8, line, " ");
            if (std.mem.startsWith(u8, cleaned, "<Hit_def>")) {
                def = try allocator.dupe(u8, cleaned[9 .. cleaned.len - 10]);
            }

            if (std.mem.startsWith(u8, cleaned, "<Hsp_identity>")) {
                size = try std.fmt.parseInt(usize, cleaned[14 .. cleaned.len - 15], 10);
            }
        }
        defer allocator.free(def); // only free if it is defined in the first place
        return BestLocalBlastHit.init(allocator, size, def) catch return null;
    }
    return null;
}

// This is code to integrate for analyzing all plasmids at once...
// const KpFolder = struct {
//     name: []const u8,
//     seqFile: ?[]const u8 = null,
// };

// pub fn main(init: std.process.Init) !void {
//     const io = init.io;
// //    const allocator = init.gpa;

//     const args = try init.minimal.args.toSlice(init.arena.allocator());
//     if (args.len != 2) {
//         std.debug.print("Usage: {s} <directory>\n", .{args[0]});
//         return;
//     }

//     const directory = args[1];
//     const cwd = std.Io.Dir.cwd();
//     const dir = cwd.openDir(io, directory, .{ .iterate = true }) catch |err| {
//         return err;
//     };
//     var it = dir.iterate();
//     while (try it.next(io)) |entry| {
//         switch (entry.kind) {
//             .file => {
//                 const start = std.mem.startsWith(u8, entry.name, "pTAN");
//                 const end = std.mem.endsWith(u8, entry.name, ".fa");
//                 const length = (entry.name.len == 10);
//                 if (start and end and length) {
//                     std.debug.print("{s}\n{s}\n", .{entry.name, @tagName(entry.kind)});
//                     }
//             },
//             .directory => {
//                 const start = std.mem.startsWith(u8, entry.name, "pTAN");
//                 if (start) {
//                     const answer = try hasExpectedFile(io, directory, entry.name);
//                     if (!answer){
//                         std.debug.print("Directory: {s}\n", .{entry.name});
//                         std.debug.print("Does not have expected file.\n", .{});
//                     }
//                 }
//             },
//             else => { continue; },
//         }
//     }
// }

// fn hasExpectedFile(io: Io, directory: []const u8, folder: []const u8) !bool {
//     var buff: [1024]u8 = undefined;
//     const path = try std.fmt.bufPrint(&buff, "{s}/{s}", .{directory,folder});
//     const subdir = std.Io.Dir.cwd().openDir(io, path, .{ .iterate = true }) catch |err| {
//         return err;
//     };
//     var it = subdir.iterate();
//     while (try it.next(io)) |entry| {
//         switch (entry.kind) {
//             .file => {
//                 const start = std.mem.startsWith(u8, entry.name, "pTAN");
//                 const end = std.mem.endsWith(u8, entry.name, ".fa");
//                 const length = (entry.name.len == 10);
//                 if (start and end and length) {
//                     return true;
//                     }
//             },
//             else => { continue; },
//         }
//     }
//     return false;
// }
